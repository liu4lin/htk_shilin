/* ----------------------------------------------------------- */
/*                                                             */
/*                          ___                                */
/*                       |_| | |_/   SPEECH                    */
/*                       | | | | \   RECOGNITION               */
/*                       =========   SOFTWARE                  */ 
/*                                                             */
/*                                                             */
/* ----------------------------------------------------------- */
/* developed at:                                               */
/*                                                             */
/*      Machine Intelligence Laboratory                        */
/*      Department of Engineering                              */
/*      University of Cambridge                                */
/*      http://mi.eng.cam.ac.uk/                               */
/*                                                             */
/* ----------------------------------------------------------- */
/*         Copyright:                                          */
/*         2002-2003  Cambridge University                     */
/*                    Engineering Department                   */
/*                                                             */
/*   Use of this software is governed by a License Agreement   */
/*    ** See the file License for the Conditions of Use  **    */
/*    **     This banner notice must not be removed      **    */
/*                                                             */
/* ----------------------------------------------------------- */
/*         File: HLVRec-outP.c OutP calculation and caching    */
/* ----------------------------------------------------------- */


static void ResetOutPCache (OutPCache *cache)
{
   int i;

   
   if (cache->nStates > 0)
      for (i = 0; i <= cache->nStates; ++i)
         cache->stateT[i] = -1000;
   
   if (cache->nMix > 0)
      for (i = 0; i <= cache->nMix; ++i)
         cache->mixT[i] = -1000;

   if (cache->nRWgts > 0)
      for (i = 0; i <= cache->nRWgts; ++i)
         cache->rwgtT[i] = -1000;

   if (cache->nTVWRs >0)
	   for (i = 0; i <= cache->nTVWRs; ++i)
		   cache->tvwrT[i] = -1000;
   cache->cacheHit = cache->cacheMiss = 0;
}

static OutPCache *CreateOutPCache (MemHeap *heap, HMMSet *hset, int block)
{
   OutPCache *cache;

   cache = New (heap, sizeof (OutPCache));

   cache->block = block;
   cache->nStates = hset->numSharedStates;/*numSharedStates is modified due to the possible inconsistent num mix per state in ConvertHSet@HLVModel.c, so it is not the traditional numSharedStates*/
   if(hset->tiedTVWR){
	   cache->nMix=hset->numMix;
	   cache->nRWgts=hset->numRegWgt;
	   cache->nTVWRs=hset->numTVWRs;/*Note that we should NOT use numSharedTVWRs, we only support the equal num mix per state case at this moment*/
   }
   else{
	   cache->nMix = hset->numSharedMix;
	   cache->nRWgts=0;
	   cache->nTVWRs=0;
   }
   /* the sIdx values are 1..numSharedStates, thus the +1 below. Same for mIdx */
   
   cache->stateOutP = cache->mixOutP = cache->lwjmts=NULL;
   if (cache->nStates > 0) {
      cache->stateT = (int *) New (heap, (cache->nStates + 1) * sizeof (int));
      cache->stateOutP = (LogFloat *) New (heap, (cache->nStates + 1) * cache->block * sizeof (LogFloat));
   }
   if (cache->nMix > 0) {
      cache->mixT = (int *) New (heap, (cache->nMix + 1) * sizeof (int));
      cache->mixOutP = (LogFloat *) New (heap, (cache->nMix + 1) * cache->block * sizeof (LogFloat));
   }
   if (cache->nRWgts>0){
	   cache->rwgtT = (int *) New (heap, (cache->nRWgts + 1) * sizeof (int));
	   cache->tvwrT = (int *) New (heap, (cache->nTVWRs + 1) * sizeof (int));
	   cache->lwjmts = (LogFloat *) New (heap, (cache->nRWgts + 1) * cache->block * sizeof (LogFloat));
   }

   return cache;
}

/* SOutP_ID_mix_Block: returns log prob of stream s of observation x */
LogFloat SOutP_ID_mix_Block(HMMSet *hset, int s, Observation *x, StreamElem *se)
{
   int vSize;
   LogDouble px;
   MixtureElem *me;
   MixPDF *mp;
   Vector v;
   LogFloat wt;

   assert (hset->hsKind == PLAINHS && hset->hsKind == SHAREDHS);
   
   v = x->fv[s];
   vSize = VectorSize(v);
   assert (vSize == hset->swidth[s]);
   me = se->spdf.cpdf+1;
   if (se->nMix == 1){     /* Single Mixture Case */
      mp = me->mpdf; 
      assert (mp->ckind == INVDIAGC);
      /*       px = IDOutP(v,vSize,mp); */
      {
         int i;
         float sum;
         float *mean, *ivar;

         mean = mp->mean;
         ivar = mp->cov.var;

         sum = mp->gConst;
         for (i=1 ; i <= vSize; i++) {
            sum += (v[i] - mean[i]) * (v[i] - mean[i]) * ivar[i];
         }
         px = -0.5*sum;
      }


      return px;
   } else {             /* Multi Mixture Case */
      LogDouble bx = LZERO;                   
      int m;

      for (m=1; m<=se->nMix; m++,me++) {
         wt = MixLogWeight(hset,me->weight);
         if (wt>LMINMIX) {  
            mp = me->mpdf; 
            /*       px = IDOutP(v,vSize,mp);   */
            {
               int i;
               float sum,xmm;
               
               sum = mp->gConst;
               for (i=1;i<=vSize;i++) {
                  xmm = v[i] - mp->mean[i];
                  sum += xmm*xmm*mp->cov.var[i];
               }
               px = -0.5*sum;
            }
            
            bx = LAdd(bx,wt+px);
         }
      }
      return bx;
   }
   return LZERO;;
}

#if 0           /* old OutPBlock()  copes with streams and non-diag outp's */
static void OutPBlock (DecoderInst *dec, Observation **obsBlock, 
                       int n, HLink hmm, int state, float acScale, LogFloat *outP)
{
   int i;
   
#if 0
   for (i = 0; i < n; ++i) {
      outP[i] = OutP (obsBlock[i], hmm, state);
   }

#else
   StateInfo *si;
   StreamElem *se;
   int s, S = obsBlock[0]->swidth[0];
   
   si = (hmm->svec+state)->info;
   se = si->pdf+1;
   
   if (S == 1 && !si->weights) {
      
      for (i = 0; i < n; ++i) {
         outP[i] = OutP_lv (dec->si, hmm->svec[state].info->sIdx, &obsBlock[i]->fv[1][1]);
#if 1   /* sanity checking */
         {
            LogFloat soutp;
            soutp = SOutP (hmm->owner, 1, obsBlock[i], se);
            assert (fabs (outP[i] - soutp) < 0.01);
         }
#endif
      }
   }
   else {       /* multi stream */
      Vector w;

      for (i = 0; i < n; ++i)
         outP[i] = 0.0;
      
      w = si->weights;
      for (s = 1; s <= S; s++, se++)
         for (i = 0; i < n; ++i)
            outP[i] += w[s] * SOutP (hmm->owner, s, obsBlock[i], se);
   }

   /* acoustic scaling */
   if (acScale != 1.0)
      for (i = 0; i < n; ++i)
         outP[i] *= acScale;
#endif
}

#endif


/* cOutP

     caching version of OutP from HModel. This only caches only on a state 
     level, not on a mixture level. 
*/
static LogFloat cOutP (DecoderInst *dec, Observation *x, HLink hmm, int state)
{
   int sIdx, n, zIdx, m, i, rIdx, n1, n2, t,s;
   LogFloat outP,bx;
   OutPCache *cache;
   StreamElem *ste;
   MixtureElem *me;
   RegWgt **rwgts;
   HMMSet* hset=dec->hset;
   HDFInfo* hdfInfo=&dec->hdfInfo;
   assert (x == dec->obsBlock[0]);
   cache = dec->outPCache;
   /*char *name=HMMPhysName(hset,hmm);*/
   sIdx = hmm->svec[state].info->sIdx;
   ste=hmm->svec[state].info->pdf+1;

   assert (sIdx >= 0);
   assert (sIdx < cache->nStates);

   if(hset->tiedTVWR){
	   zIdx = hmm->zvec[state].info->zIdx;
	   assert (zIdx >= 0);
	   assert (zIdx < cache->nTVWRs);
	   n1=dec->frame-cache->tvwrT[zIdx];/*this is used to track cached TVWR, a collection of regression weights*/
	   assert (n1>=0);
	   for(s=1;s<=hset->swidth[0];s++){
		   UpdateRegWgts(ste, hmm->zvec[state].info, s);
		   if(n1>=cache->block){
			   rwgts=hmm->zvec[state].info->tstr[s].rwgts;
			   for(m=1;m<=ste->nMix;m++){
				   rIdx=rwgts[m]->rIdx;
				   assert (rIdx >= 0);
				   assert (rIdx < cache->nRWgts);
				   n2=dec->frame-cache->rwgtT[rIdx];
				   assert (n2>=0);
				   if(n2>=cache->block){
					   for(i=0;i<cache->block;i++){
						   t=dec->frame+i-hdfInfo->StartTime;
						   /*MulDotProduct(dec->hdfInfo.ctxPos[t], rwgts[m]->regparm, hset->regwsize, rwgts[m]->wjmt);*/
						   switch(hset->rkind){
						   case NONE_CTX_REG:
						   case TEMP_CTX_REG:
							   TempCtxExpReg(hdfInfo->frmPos[s], hdfInfo->ctxPos[s], hdfInfo->ctxLoc[s], rwgts[m]->regparm, hset->pwidth[s], rwgts[m]->wjmt, t, hset->temp_ctx_size, hset->priors?hset->priors[s]:NULL);
							   break;
						   case SPAT_CTX_REG:
							   SpatCtxExpReg(hdfInfo->frmPos, rwgts[m]->regparm, hset->pwidth, rwgts[m]->wjmt, t, hset->priors);
							   break;
						   case SPTP_CTX_REG:
							   SptpCtxExpReg(hdfInfo->frmPos, hdfInfo->ctxPos, hdfInfo->ctxLoc, rwgts[m]->regparm, hset->pwidth, rwgts[m]->wjmt, t, hset->temp_ctx_size, hset->priors);
							   break;
						   case HIGH_SPA_REG:
						   case TEMP_SPA_REG:
							   TempCtxSpaReg(hdfInfo->frmPos[s], hdfInfo->ctxPos[s], hdfInfo->ctxLoc[s], rwgts[m]->regparm, rwgts[m]->wjmt, t, hset->temp_ctx_size, hset->confuTbl[rwgts[m]->spIdx], hset->priors?hset->priors[s]:NULL);
							   break;
						   }
						   cache->lwjmts[rIdx*cache->block+i]=rwgts[m]->wjmt[0];/*last bug found, miss taking log*/
						   rwgts[m]->cached_t=t;
					   }
					   cache->rwgtT[rIdx]=dec->frame;
				   }
			   }
		   }
		   cache->tvwrT[zIdx]=dec->frame;
	   }
   }

   n = dec->frame - cache->stateT[sIdx];

   assert (n >= 0);

   if (n < cache->block) {
      outP = cache->stateOutP[sIdx * cache->block + n];
      if(hset->tiedTVWR){/*re-calculate stateOutP from cached mixOutP and wjmt, wjmt may need to re-calculate from the data*/
    	  outP=0;
    	  me=ste->spdf.cpdf+1;
    	  for(s=1;s<=hset->swidth[0];s++){
        	  rwgts=hmm->zvec[state].info->tstr[s].rwgts;
        	  bx=LZERO;
        	  for(m=1;m<=ste->nMix;m++,me++){
        		  bx=LAdd(bx,cache->lwjmts[rwgts[m]->rIdx*cache->block+n]+cache->mixOutP[me->mpdf->mIdx*cache->block+n]);/*duplicate standard outp for multiple posterior stream*/
        	  }
        	  outP+=bx;
    	  }
    	  outP *=dec->acScale;
      }
      ++cache->cacheHit;
#if 0
      /* the following is *very* expensive, it effectively disables the cache,
         use only for sanity checking! */
      assert (outP == dec->acScale * OutP (x, hmm, state));
#endif
   }
   else {
      ++cache->cacheMiss;
      if (!cache->mixOutP||hset->tiedTVWR) {     /* don't bother caching mixtures */
         /* #### handle boundary case where we don't have cache->block obs left */

         if (!dec->si->useHModel) 
            OutPBlock (dec->si, &dec->obsBlock[0], cache->block,
                       sIdx, dec->acScale, &cache->stateOutP[sIdx * cache->block], cache->mixOutP,
                       dec->frame);
         else{
            OutPBlock_HMod (dec->si, &dec->obsBlock[0], cache->block,
                            sIdx, dec->acScale, &cache->stateOutP[sIdx * cache->block], cache->mixOutP,
                            dec->frame);
         }
            
         cache->stateT[sIdx] = dec->frame;
         outP = cache->stateOutP[sIdx * cache->block];

#if 0   /* sanity checking for OutPBlock */
         {
            LogFloat safe_outP;
            safe_outP = dec->acScale * OutP (x, hmm, state);
            assert (fabs (outP - safe_outP) < 0.01);
         }
#endif
      }
      else {            /* cache mixtures (e.g. for soft-tied systems) */
         abort ();
         /*
x      outP = OutP (x, hmm, state);
x      dec->cacheOutP[sIdx] = outP;
x      CACHE_FLAG_SET(dec, sIdx);
         */
      }
   }      

   return outP;
}



/* outP caclulation for USEHMODEL=T case  */


/*******************************************************************************/
/*  outP calculation from HModel.c and extended for new adapt code */


static LogFloat SOutP_HMod (HMMSet *hset, int s, Observation *x, StreamElem *se, LogFloat* mixOutP, int id, int block, int block_i)
{
   int m;
   LogFloat bx,px,wt,det=0,swt;
   MixtureElem *me;
   MixPDF *mp;
   Vector v,otvs;
   
   /* Note hset->kind == SHAREDHS */
   assert (hset->hsKind == SHAREDHS);

   v=x->fv[s];
  if(hset->fMPE){
	  v=ApplyFprojection((HDFInfo*)hset->hdfInfo,id);
  }
  if(hset->TVWR){
	  ApplyWprojection(se, (HDFInfo*)hset->hdfInfo, id, s);
  }
   me=se->spdf.cpdf+1;
   if (se->nMix==1){     /* Single Mixture Case */
	  px= hset->disableGMM?0:MOutP(ApplyCompFXForm(me->mpdf,v,inXForm,&det,id),me->mpdf);
      px += det;
      if(hset->tiedTVWR){
    	  mixOutP[me->mpdf->mIdx*block+block_i]=px;
      }
      swt=0;
      bx = MixLogWeight(hset,me->weight)+px;
   } else if (!pde) {
      bx=LZERO;                   /* Multi Mixture Case */
      float count=0;
      for (m=1; m<=se->nMix; m++,me++) {
    	 swt = MixLogWeight(hset,me->sweight);
    	 wt = hset->TVWR?MixLogWeight(hset,me->weight):swt;
         if (swt>LMINMIX){
        	 count++;
        	 px= hset->disableGMM?0:MOutP(ApplyCompFXForm(me->mpdf,v,inXForm,&det,id),me->mpdf);
        	 px += det;
        	 bx=LAdd(bx,wt+px);
         }
         if(mixOutP){
        	 mixOutP[me->mpdf->mIdx*block+block_i]=swt>LMINMIX?swt+px:LZERO;
         }
      }
   } else {   /* Partial distance elimination */
      wt = MixLogWeight(hset,me->weight);
      mp = me->mpdf;
      otvs = ApplyCompFXForm(mp,v,inXForm,&det,id);
      px = IDOutP(otvs,VectorSize(otvs),mp);
      bx = wt+px+det;
      for (m=2,me=se->spdf.cpdf+2;m<=se->nMix;m++,me++) {
          wt = MixLogWeight(hset,me->weight);
          swt = MixLogWeight(hset,me->sweight);
          if (swt>LMINMIX){
        	 mp = me->mpdf;
        	 otvs = ApplyCompFXForm(mp,v,inXForm,&det,id);
        	 if (PDEMOutP(otvs,mp,&px,bx-wt-det) == TRUE){
        		 px+=det;
        		 bx = LAdd(bx,wt+px);
        	 }
         }
      }
   }
   return bx;
}

LogFloat POutP_HModel (HMMSet *hset,Observation *x, StateInfo *si, LogFloat* mixOutP, int id, int block, int block_i)
{
   LogFloat bx;
   StreamElem *se;
   Vector w;
   int s,S = x->swidth[0];
   
   if (S==1 && si->weights==NULL)
      return SOutP_HMod(hset,1,x,si->pdf+1, mixOutP, id, block, block_i);
   bx=0.0; se=si->pdf+1; w = si->weights;
   for (s=1;s<=S;s++,se++){
      bx += w[s]*SOutP_HMod(hset,s,x,se, mixOutP, id, block, block_i);
   }
   return bx;
}

void OutPBlock_HMod (StateInfo_lv *si, Observation **obsBlock, 
                int n, int sIdx, float acScale, LogFloat *outP, LogFloat *mixOutP, int id)
{
   int i;

   assert  (si->useHModel);
   
   for (i = 0; i < n; ++i) {
      outP[i] = POutP_HModel (si->hset, obsBlock[i], si->si[sIdx], mixOutP, id+i, n, i);
   }
   /* acoustic scaling */
   if (acScale != 1.0)
      for (i = 0; i < n; ++i)
         outP[i] *= acScale;
}
