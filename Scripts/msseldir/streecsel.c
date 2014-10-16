/**********  segtre_mig.c **********************************
*
*	This subroutine uses a Monte Carlo algorithm described in
*	Hudson,R. 1983. Theor.Pop.Biol. 23: 183-201, to produce
*	a history of a random sample of gametes under a neutral
*	Wright-Fisher model with recombination and geographic structure. 
*	Input parameters
*	are the sample size (nsam), the number of sites between
*	which recombination can occur (nsites), and the recombination
*	rate between the ends of the gametes (r). The function returns
*	nsegs, the number of segments the gametes were broken into
*	in tracing back the history of the gametes.  The histories of
*	these segments are passed back to the calling function in the
*	array of structures seglst[]. An element of this array,  seglst[i],
* 	consists of three parts: (1) beg, the starting point of
*	of segment i, (2) ptree, which points to the first node of the
*	tree representing the history of the segment, (3) next, which
*	is the index number of the next segment.
*	     A tree is a contiguous set of 2*nsam nodes. The first nsam
*	nodes are the tips of the tree, the sampled gametes.  The other
*	nodes are the nodes ancestral to the sampled gametes. Each node
*	consists of an "abv" which is the number of the node ancestral to
*	that node, an "ndes", which is not used or assigned to in this routine,
*	and a "time", which is the time (in units of 4N generations) of the
*	node. For the tips, time equals zero.
*	Returns a pointer to an array of segments, seglst.

**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mssel.h"
#define NL putchar('\n')
#define size_t unsigned

#define MIN(x, y) ( (x)<(y) ? (x) : (y) )
#define MAX(x, y) ( (x)>(y) ? (x) : (y) )

#define ERROR(message) fprintf(stderr,message),NL,exit(1)

#define SEGINC 80 

extern int flag;

int nchrom, begs, nsegs, selspot;
long nlinks ;
static int *nnodes = NULL ;  
double t, cleft , pc, lnpc, *freqder  ;

static unsigned seglimit = SEGINC ;
static unsigned maxchr ;
static int *config[2] = {NULL}, **inconfig ;

struct seg{
	int beg;
	int end;
	int desc;
	};

struct chromo{
	int nseg;
	int pop;
        int allele;
	struct seg  *pseg;
	} ;

static struct chromo *chrom = NULL ;

struct node{
	int abv;
	int ndes;
	float time;
	} *ptree1, *ptree2;

struct segl {
	int beg;
	struct node *ptree;
	int next;
	}  ;
static struct segl *seglst = NULL ;

	struct segl *
segtre_mig(struct c_params *cp, int *pnsegs ) 
{
	int i, j, k,  dec, pop, pop2, c1, c2, ind, rchrom  ;
	int migrant, source_pop, flagint ;
	double  ran1(), sum, x, ttemp, rft, clefta,  tmin, p  ;
	double prec, cin,  prect,  mig, ran, coal_prob, prob, rdum , arg ;
	char  event ;
	int re(), xover(), cinr(), cleftr(), eflag, cpop, ic, gap, gapleft, gapright  ;
	int nsam, npop, nsites, allele, callele, freqcount ,pop1, ipop, mallele  ;
	double r,  f, rf,  track_len,  **migm ;
	double *size, *alphag, *tlast,  rbet, fr, tnextfreq  ;
	double coal_proba, coal_probd, mig2, sumfreq ;
	struct devent *nextevent ;
	struct seg *pseg ;
	void pick2_chrom( ) ;

	nsam = cp->nsam;
	npop = cp->npop;
	nsites = cp->nsites;
	inconfig = cp->config;
	r = cp->r ;
	f = cp->f ;
	track_len = cp->track_len ;
	migm = (double **)malloc( (unsigned)npop*sizeof(double *) ) ;
	for( i=0; i<npop; i++) {
	  migm[i] = (double *)malloc( (unsigned)npop*sizeof( double) ) ;
	  for( j=0; j<npop; j++) migm[i][j] = (cp->mig_mat)[i][j] ;
	  }
	nextevent = cp->deventlist ;
	selspot = cp->selspot ;
	
/* Initialization */
	if( chrom == NULL ) {
	   maxchr = nsam + 20 ;
	   chrom = (struct chromo *)malloc( (unsigned)( maxchr*sizeof( struct chromo) )) ;
	  if( chrom == NULL ) perror( "malloc error. segtre");
	  }
	if( nnodes == NULL ){
		nnodes = (int*) malloc((unsigned)(seglimit*sizeof(int)))  ;
		if( nnodes == NULL ) perror("malloc error. segtre_mig");
		}
	if( seglst == NULL ) {
		seglst = (struct segl *)malloc((unsigned)(seglimit*sizeof(struct segl)) ) ;
		if( seglst == NULL ) perror("malloc error. segtre_mig.c 2");
		}

	if( config[0] == NULL ){
	    config[0] = (int *)malloc( (unsigned)((npop+1)*sizeof(int) )) ;
	    config[1] = (int *)malloc( (unsigned)((npop+1)*sizeof(int) )) ;
	}
	if( config[0] == NULL ) perror("malloc error. segtre.");
	size = (double *)malloc( (unsigned) ((npop)*sizeof(double) )) ;
	if( size == NULL ) perror("malloc error. segtre.");
	alphag = (double *)malloc( (unsigned) ((npop)*sizeof(double) )) ;
	if( alphag == NULL ) perror("malloc error. segtre.");
	tlast = (double *)malloc( (unsigned) ((npop)*sizeof(double) )) ;
	if( tlast == NULL ) perror("malloc error. segtre.");
	for(pop=ind=0;pop<npop;pop++){
	   size[pop] = (cp->size)[pop] ;
	   alphag[pop] = (cp->alphag)[pop] ;
	   tlast[pop] = 0.0 ;
	   for(allele = 0; allele < 2 ; allele++){
	        config[allele][pop] = inconfig[allele][pop];
		    for(j=0; j<inconfig[allele][pop];j++,ind++) {
			
			  chrom[ind].nseg = 1;
			  if( !(chrom[ind].pseg = (struct seg*)malloc((unsigned)sizeof(struct seg)) ))
			      ERROR("calloc error. se1");

			   (chrom[ind].pseg)->beg = 0;
			   (chrom[ind].pseg)->end = nsites-1;
			   (chrom[ind].pseg)->desc = ind ;
			   chrom[ind].pop = pop ;
			   chrom[ind].allele = allele ;
			}
	   }
	}
	seglst[0].beg = 0;
	if( !(seglst[0].ptree = (struct node *)calloc((unsigned)(2*nsam),sizeof(struct node)) ))
		 perror("calloc error. se2");



	nnodes[0] = nsam - 1 ;
	nchrom=nsam;
	nlinks = ((long)(nsam))*(nsites-1) ;
	nsegs=1;
	t = 0.;
	r /= (nsites-1);
	if( f > 0.0 ) 	pc = (track_len -1.0)/track_len ;
	else pc = 1.0 ;
	lnpc = log( pc ) ;
	cleft = nsam* ( 1.0 - pow( pc, (double)(nsites-1) ) ) ;
	if( r > 0.0 ) rf = r*f ;
	else rf = f /(nsites-1) ;
	rft = rf*track_len ;
	flagint = 0 ;
	freqcount = 1 ;
	tnextfreq = cp->tfreqder[freqcount] ;
	freqder = cp->infreqder[0] ;

/* Main loop */

	while( nchrom > 1 ) {
		prec = nlinks*r;
		cin = nlinks*rf ;
		clefta = cleft*rft ;
		prect = prec + cin + clefta ;
		rbet = 0 ;
		for( ic=0; ic<nchrom ; ic++) {  
		  pseg = chrom[ic].pseg;
		  gapleft =  pseg->beg - selspot ;
		  gapright =  selspot - ((pseg+chrom[ic].nseg -1)->end ) ;
		  gap = MAX( gapleft, gapright ) ;
		  if( gap > 0 ){
		    if( chrom[ic].allele == 0 ) fr =  freqder[ chrom[ic].pop ] ;
		    else fr = 1.0 - freqder[ chrom[ic].pop ] ;
		    rbet += fr*(gap*r + rft*(1.0-pow(pc,(double)gap))) ;
		  }
		}
		sumfreq = 0.0 ;
		for( pop =0; pop < cp->maxpop; pop++) sumfreq += freqder[pop] ;

		eflag = 0;
		tmin = 99999.0 ;
	/*	fprintf(stderr,"%lf %d %d %lf\n", t,config[0][0],config[1][0],freqder[0]); */

		mig = 0.0;  /* calc mig rate, if any freqs are zero, instaneous migration possible.  */
		mig2 = 0.0 ;
		for( pop1 =0; pop1 <npop; pop1++ ) {
		    if( config[0][pop1] > 0 ){
		      if ( freqder[pop1] < 1.0 ){
		        for( pop2 = 0; pop2<npop ; pop2++){
		          if( pop2 != pop1 )
			    mig += config[0][pop1]*migm[pop1][pop2]*(1.-freqder[pop2])/(1.-freqder[pop1]) ;
			    }
		      }
		      else if( (freqder[pop1] == 1) && (config[0][pop1] == 1) ){
			    mig2 = 0. ;
			    for( pop2 = 0; pop2 < npop; pop2++)
			       if( pop2 != pop1 ) mig2 +=  config[0][pop1]*migm[pop1][pop2]*(1.-freqder[pop2]);
			    if( mig2 > 0. ){
			       eflag = 1 ;
			       event = 'i';
			       tmin = 0.0 ;
			       ipop = pop1 ;
			       mallele = 0 ;
			       mig = 0.0 ;
			       break;
			    }
		      }
		    }
		    if( config[1][pop1] > 0 ) {
		      if ( freqder[pop1] > 0.0 ){
		         for( pop2 = 0; pop2<npop ; pop2++){
		            if( pop2 != pop1 )
			      mig += config[1][pop1]*migm[pop1][pop2]*freqder[pop2]/freqder[pop1] ;
			     }
		      }
		      else if( (freqder[pop1] == 0 ) && ( config[1][pop1] == 1)){
			    if( sumfreq == 0.0 ) {
			       eflag = 1;
			       tmin = 0.0 ;
			       cpop = pop1 ;
			       event = 'M' ;
			       break;
			    }
			    mig2 = 0.0 ;
			    for( pop2 = 0; pop2 < npop; pop2++)
			    if( pop1 != pop2 ) mig2 += config[1][pop1]*migm[pop1][pop2]*freqder[pop2] ;
			    if( mig2 > 0.0 ) {
			       eflag = 1 ;
			       event = 'i';
			       tmin = 0.0 ;
			       ipop = pop1 ;
			       mallele = 1 ;
			       mig = 0.0 ;
			       break;
				}
		      }
		    }
		}

		if( (npop > 1) && ( mig == 0.0) && (mig2 == 0.0) && ( nextevent == NULL)) {
		   i = 0;
		   for( j=0; j<npop; j++) 
		     if( (config[0][j]+config[1][j]) > 0 ) i++;
		   if( i > 1 ) {
			fprintf(stderr," Infinite coalescent time. No migration.\n");
			exit(1);
		   }
		}

		if( prect > 0.0 ) {      /* cross-over or gene conversion */
		  while( (rdum = ran1() )  == 0.0 ) ;
		  ttemp = -log( rdum)/prect ;
		  if( (eflag == 0) || (ttemp < tmin ) ){
		    tmin = ttemp;
		    event = 'r' ;
		    eflag = 1;
		  }
	        }
		
		if( rbet > 0.0 ){
		  while( (rdum = ran1() )  == 0.0 ) ;
		  ttemp = -log( rdum)/rbet ;
		  if( (eflag == 0) || (ttemp < tmin ) ){
		    tmin = ttemp;
		    event = 'b' ;
		    eflag = 1;
		  }

		}

		if(mig > 0.0 ) {         /* migration   */
		  while( (rdum = ran1() ) == 0.0 ) ;
		  ttemp = -log( rdum)/mig ;
		  if( (eflag == 0) || (ttemp < tmin ) ){
		    tmin = ttemp;
		    event = 'm' ;
		    eflag = 1 ;
		  }
	        }
		
	    for(pop=0; pop<npop ; pop++) {     /* coalescent */
	       if( (config[0][pop] > 1) && ( freqder[pop] == 1 )) {
		      tmin = 0.0 ;
		      eflag = 1 ;
		      cpop = pop ;
		      callele = 0;
		      event = 'c' ;
		      break;
		   }
	       if( freqder[pop] < 1.0 )
	          coal_proba = ((double)config[0][pop])*(config[0][pop]-1.)/(1.-freqder[pop]) ;
	       else coal_proba = 0.0 ;
	       if( (config[1][pop] > 1) && ( freqder[pop] == 0 )) {
		      tmin = 0.0 ;
		      eflag = 1 ;
		      cpop = pop ;
		      callele = 1;
		      event = 'c' ;
		      break;
		   }
	       if( freqder[pop] > 0.0)
		     coal_probd = ((double)config[1][pop])*(config[1][pop]-1.)/freqder[pop] ;
	       else coal_probd = 0.0 ;
	       coal_prob = coal_proba + coal_probd ;
	       if( coal_prob > 0.0 ) {
		   while( ( rdum = ran1() )  == .0 ) ;
	  	   if( alphag[pop] == 0 ){
			 ttemp = -log( rdum )*size[pop] /coal_prob ;
		         if( (eflag == 0) || (ttemp < tmin ) ){
		            tmin = ttemp;
		            event = 'c' ;
		            eflag = 1 ;
			    cpop = pop;
			    if( ran1() < coal_proba/coal_prob ) callele = 0;
			    else callele = 1 ;
			 }
		    }
	  	   else {
		     arg  = 1. - alphag[pop]*size[pop]*exp(-alphag[pop]*(t - tlast[pop] ) )* log(rdum) / coal_prob     ;
		     if( arg > 0.0 ) {                          /*if arg <= 0,  no coalescent within interval */ 
		         ttemp = log( arg ) / alphag[pop]  ;
		         if( (eflag == 0) || (ttemp < tmin ) ){
		            tmin = ttemp;
		            event = 'c' ;
		            eflag = 1 ;
			    cpop = pop ;
			    if( ran1() < coal_proba/coal_prob ) callele = 0;
			    else callele = 1 ;
			 }
		     }

		   }
	         }		
 	      }
	    if( freqcount < cp->nfreqder ) 
		         if( (eflag == 0) || ( tnextfreq - t < tmin ) ){
		            tmin = tnextfreq-t ;
		            event = 'f' ;
		            eflag = 1 ;
	
			 }

	    if( (eflag == 0) && ( nextevent == NULL) ) {
	      fprintf(stderr,
               " infinite time to next event. Negative growth rate in last time interval or non-communicating subpops.\n");
			   fprintf(stderr," %lf %lf %lf %d %d %lf\n", tnextfreq, t, tmin, freqcount , cp->nfreqder, freqder[0] );
	      exit( 0);
	    }
	if( ( ( eflag == 0) && (nextevent != NULL))|| ( (nextevent != NULL) &&  ( (t+tmin) >=  nextevent->time)) ) {
	    t = nextevent->time ;
	    switch(  nextevent->detype ) {
		case 'N' :
		   for(pop =0; pop <npop; pop++){
			 size[pop]= nextevent->paramv ;
			 alphag[pop] = 0.0 ;
			}
		   nextevent = nextevent->nextde ;
		   break;
		case 'n' :
		   size[nextevent->popi]= nextevent->paramv ;
		   alphag[nextevent->popi] = 0.0 ;
		   nextevent = nextevent->nextde ;
		   break;
		case 'G' :
		   for(pop =0; pop <npop; pop++){
		     size[pop] = size[pop]*exp( -alphag[pop]*(t - tlast[pop]) ) ;
		     alphag[pop]= nextevent->paramv ;
		     tlast[pop] = t ;
		   }
		   nextevent = nextevent->nextde ;
		   break;
		case 'g' :
		     pop = nextevent->popi ;
		     size[pop] = size[pop]*exp( - alphag[pop]*(t-tlast[pop]) ) ;
		     alphag[pop]= nextevent->paramv ;
		     tlast[pop] = t ;
		     nextevent = nextevent->nextde ;
		     break;
		case 'M' :
		   for(pop =0; pop <npop; pop++)
		     for( pop2 = 0; pop2 <npop; pop2++) migm[pop][pop2] = (nextevent->paramv) /(npop-1.0) ;
		   for( pop = 0; pop <npop; pop++)
		     migm[pop][pop]= nextevent->paramv ;
		   nextevent = nextevent->nextde ;
		   break;
		case 'a' :
		   for(pop =0; pop <npop; pop++)
		     for( pop2 = 0; pop2 <npop; pop2++) migm[pop][pop2] = (nextevent->mat)[pop][pop2]  ;
		   nextevent = nextevent->nextde ;
		   break;
		case 'm' :
		  i = nextevent->popi ;
		  j = nextevent->popj ;
		  migm[i][i] += nextevent->paramv - migm[i][j];
		   migm[i][j]= nextevent->paramv ;
		   nextevent = nextevent->nextde ;
		   break;
	        case 'j' :         /* merge pop i into pop j  (join) */
		  i = nextevent->popi ;
		  j = nextevent->popj ;
		  config[0][j] += config[0][i] ;
		  config[1][j] += config[1][i] ;
		  config[0][i] = config[1][i] = 0 ;
		  for( ic = 0; ic<nchrom; ic++) if( chrom[ic].pop == i ) chrom[ic].pop = j ;
		/*  the following was added 19 May 2007 */
		  for( k=0; k < npop; k++){
		     if( k != i) {
		        migm[k][k] -= migm[k][i] ;
		        migm[k][i] = 0. ;
			 }
		   }
		/* end addition */
		   nextevent = nextevent->nextde ;
		   break;
	        case 's' :         /*split  pop i into two;p is the proportion from pop i, and 1-p from pop n+1  */
		  i = nextevent->popi ;
		  p = nextevent->paramv ;

		  npop++;
		  config[0] = (int *)realloc( config[0], (unsigned)(npop*sizeof( int) )); 
		  config[1] = (int *)realloc( config[1], (unsigned)(npop*sizeof( int) )); 
		  size = (double *)realloc(size, (unsigned)(npop*sizeof(double) ));
		  alphag = (double *)realloc(alphag, (unsigned)(npop*sizeof(double) ));
		  tlast = (double *)realloc(tlast,(unsigned)(npop*sizeof(double) ) ) ;
		  tlast[npop-1] = t ;
		  size[npop-1] = 1.0 ;
		  alphag[npop-1] = 0.0 ;
		  migm = (double **)realloc(migm, (unsigned)(npop*sizeof( double *)));
		  for( j=0; j< npop-1; j++)
			 migm[j] = (double *)realloc(migm[j],(unsigned)(npop*sizeof(double)));
		  migm[npop-1] = (double *)malloc( (unsigned)(npop*sizeof( double) ) ) ;
		  for( j=0; j<npop; j++) migm[npop-1][j] = migm[j][npop-1] = 0.0 ;
		  config[0][npop-1] = config[1][npop-1] = 0 ;
		  config[0][i] = config[1][i] = 0 ;
		  for( ic = 0; ic<nchrom; ic++){
		    if( chrom[ic].pop == i ) {
		      allele = chrom[ic].allele ;
		      if( allele == 0 ){
			if( (1.0 - freqder[i]) > 0 )
			  prob = (1-p)*(1.0- freqder[npop-1])/(1.0-freqder[i] );
			else prob = 1.0 ;
		      }
		      else {
			if( freqder[i] > 0 )
			  prob = (1-p)*freqder[npop-1]/freqder[i] ;
			else prob = 1.0 ;
		      }
		      if( ran1() > prob ) config[ allele ][i]++; 
		      else {
			 chrom[ic].pop = npop-1 ;
			 config[ allele][npop-1]++;
		      }
		    }
		  }
		  freqder[i] = ( freqder[i] - (1.0-p)*freqder[npop-1] )/p ;
		   nextevent = nextevent->nextde ;
		   break;
		}
 	   } 
	else {
/*	fprintf(stderr,"t: %lf tmin: %lf %c  %d %d %d %d %d %d  mig: %lf  mig2: %lf\n", t,tmin,event, config[0][0], config[1][0],config[0][1], 
	   config[1][1],config[0][2], config[1][2], mig, mig2);  */
/*	   	for( ic = 0; ic< nchrom; ic++ ) fprintf(stderr," pop: %d  allele: %d\n", chrom[ic].pop, chrom[ic].allele );
	fprintf(stderr,"\n\n"); */

		   t += tmin ;	
		   if( event == 'r' ) {   
		      if( (ran = ran1()) < ( prec / prect ) ){ /*recombination*/
		     	  rchrom = re(nsam);
		      }
		      else if( ran < (prec + clefta)/(prect) ){    /*  cleft event */
			 rchrom = cleftr(nsam);
		      }
		      else  {         /* cin event */
			 rchrom = cinr(nsam,nsites);
		      }
		   }

		   else if ( event == 'b' ){
		     x = rbet*ran1();
		     sum = 0.0 ;
		     for( ic=0; ic<nchrom ; ic++) {
		        pseg = chrom[ic].pseg;
		        gapleft =  pseg->beg - selspot ;
		        gapright =  selspot - ((pseg+chrom[ic].nseg -1)->end ) ;
		        gap = MAX( gapleft, gapright ) ;
		        if( gap > 0 ){
		           if( chrom[ic].allele == 0 ) fr =  freqder[ chrom[ic].pop ] ;
		           else fr = 1.0 - freqder[ chrom[ic].pop ] ;
		           sum += fr*(gap*r + rft*(1.0-pow(pc,(double)gap))) ;
			   if( x < sum ) break;
		         }
		     }
		     config[ chrom[ic].allele][chrom[ic].pop] -= 1 ;
		     config[ 1-chrom[ic].allele][chrom[ic].pop] += 1 ;
		     chrom[ic].allele = 1 - chrom[ic].allele ;
		   }

	           else if ( event == 'm' ) {  /* migration event */
			x = mig*ran1();
			sum = 0.0 ;
			for( ind =0; ind < nchrom; ind++){
			  ipop = chrom[ind].pop ;
			  allele = chrom[ind].allele ; /* in the following need to check denominator zero cases */
			  for( pop = 0; pop < npop; pop++){
			    if( ipop != pop){
			      if( ( allele == 0) && ( freqder[ipop] < 1.0 ) ) sum += migm[ipop][pop]*(1.-freqder[pop])/(1.-freqder[ipop]) ;
			      else if( freqder[ipop] > 0. )  sum += migm[ipop][pop]*freqder[pop]/freqder[ipop] ;
			      if( x < sum ) goto gotmigrant ;
			    }
			  }
			}
		   gotmigrant:
			migrant = ind ;
			source_pop = pop ;
			  config[allele][chrom[migrant].pop] -= 1;
			  config[allele][source_pop] += 1;
			  chrom[migrant].pop = source_pop ;
	           }

		   else if( event == 'f' ){  /* change in selected allele frequency */
		     freqder = cp->infreqder[freqcount];
			 freqcount++ ;
			 if( freqcount < cp->nfreqder)  tnextfreq = cp->tfreqder[freqcount] ;
			 else tnextfreq = 99999. ;
		   }
		   else if( event == 'M' ){
		     config[1][cpop]--;
		     config[0][cpop]++;
		     for( ic=0; ic<nchrom; ic++){
		       if( (chrom[ic].pop== cpop)&&( chrom[ic].allele == 1) ) {
			 chrom[ic].allele = 0; 
			 break;
		       }
		     }
		   }
		   else if( event == 'i'){
		     x = mig2*ran1();
		     sum = 0.0;
		     for( pop =0; pop < npop; pop++){
		       if( mallele == 1 ) fr = freqder[pop];
		       else fr = 1.0 - freqder[pop] ;
		       if( pop != ipop ) sum += config[mallele][ipop]*migm[ipop][pop]*fr ; 
		       if ( x < sum ) break;
		     }
		     config[mallele][ipop]--;
		     config[mallele][pop]++;
		     for( ic = 0; ic < nchrom; ic++)
		       if( (chrom[ic].pop == ipop) && ( chrom[ic].allele == mallele ) ){
		           chrom[ic].pop = pop ;
		           break;
		       }
		     
		   }
		   else { 								 /* coalescent event */
			/* pick the two, c1, c2  */
		       pick2_chrom( cpop, config[callele],callele, &c1,&c2);  /* c1 and c2 are chrom's to coalesce */
			dec = ca(nsam,nsites,c1,c2 );
			config[callele][cpop] -= dec ;
		   }
		 }
	     }  
	*pnsegs = nsegs ;
	free( size ) ;
	free( alphag );
	free( tlast );
	for( i=0; i<npop; i++) free ( migm[i] ) ;
	free( migm ) ;
	return( seglst );
}

/******  recombination subroutine ***************************
   Picks a chromosome and splits it in two parts. If the x-over point
   is in a new spot, a new segment is added to seglst and a tree set up
   for it.   ****/


	int
re(nsam)
	int nsam;
{
	struct seg *pseg ;
	int  el, lsg, lsgm1,  ic,  is, spot, newallele, tchrom;
	double ran1();


/* First generate a random x-over spot, then locate it as to chrom and seg. */

	spot = nlinks*ran1() + 1.;

    /* get chromosome # (ic)  */

	for( ic=0; ic<nchrom ; ic++) {
		lsg = chrom[ic].nseg ;
		lsgm1 = lsg - 1;
		pseg = chrom[ic].pseg;
		el = ( (pseg+lsgm1)->end ) - (pseg->beg);
		if( spot <= el ) break;
		spot -= el ;
		}
	is = pseg->beg + spot -1;
	xover(nsam, ic, is);
	newallele  =  (ran1() < freqder[chrom[ic].pop] ) ? 1 : 0  ;
	tchrom =  (is < selspot) ? ic : nchrom-1  ;
	chrom[nchrom-1].allele = chrom[ic].allele ;
	chrom[tchrom].allele = newallele ;
	config[ newallele][ chrom[ic].pop] += 1 ;
  
	return(ic);	
}

	int
cleftr( int nsam)
{
	struct seg *pseg ;
	int  ic,  is,  newallele, tchrom ;
	double ran1(), x, sum, len  ;
	int xover(int, int, int);

	while( (x = cleft*ran1() )== 0.0 )  ;
	sum = 0.0 ;
	ic = -1 ;
	while ( sum < x ) {
		sum +=  1.0 - pow( pc, links(++ic) )  ;
		}
	pseg = chrom[ic].pseg;
	len = links(ic) ;
	is = pseg->beg + floor( 1.0 + log( 1.0 - (1.0- pow( pc, len))*ran1() )/lnpc  ) -1  ;
	xover( nsam, ic, is);
	newallele  =  (ran1() < freqder[chrom[ic].pop]) ? 1 : 0  ;
	tchrom =  (is < selspot) ? ic : nchrom-1  ; 
	chrom[nchrom-1].allele = chrom[ic].allele ;
	chrom[tchrom].allele = newallele ;
	config[ newallele][ chrom[ic].pop] += 1 ;
	return( ic) ;
}

	int
cinr( int nsam, int nsites)
{
	struct seg *pseg ;
	int len,  el, lsg, lsgm1,  ic,  is, spot, endic, newallele,tchrom ;
	double ran1();
	int xover(), ca() ;


/* First generate a random x-over spot, then locate it as to chrom and seg. */

	spot = nlinks*ran1() + 1.;

    /* get chromosome # (ic)  */

	for( ic=0; ic<nchrom ; ic++) {
		lsg = chrom[ic].nseg ;
		lsgm1 = lsg - 1;
		pseg = chrom[ic].pseg;
		el = ( (pseg+lsgm1)->end ) - (pseg->beg);
		if( spot <= el ) break;
		spot -= el ;
		}
	is = pseg->beg + spot -1;
	endic = (pseg+lsgm1)->end ;
	xover(nsam, ic, is);

	len = floor( 1.0 + log( ran1() )/lnpc ) ;
	if( is+len >= endic ){
	  newallele  =  (ran1() < freqder[chrom[ic].pop] ) ? 1 : 0  ;
	  if( (is < selspot) && ( (is+len) >= selspot ) ) tchrom = ic ;
	  else tchrom = nchrom -1  ;
	  chrom[nchrom-1].allele = chrom[ic].allele ;
	  chrom[tchrom].allele = newallele ;
	  config[ newallele][ chrom[ic].pop] += 1 ;
	  return(ic) ; 
	}
	if( is+len < (chrom[nchrom-1].pseg)->beg ){
	   ca( nsam, nsites, ic, nchrom-1);
	   if( (is < selspot) && ( (is+len) >= selspot ) ){
	     newallele  =  (ran1() < freqder[chrom[ic].pop]) ? 1 : 0  ;
	     if( newallele != chrom[ic].allele ){
	       config[ chrom[ic].allele][chrom[ic].pop] -= 1 ;
	       config[ newallele ][chrom[ic].pop] += 1 ;
	       chrom[ic].allele = newallele;
	     }   
	   }
	    return(-1) ;
	}
	xover( nsam, nchrom-1, is+len ) ;
	chrom[nchrom-1].allele = chrom[ic].allele ;
	newallele  =  (ran1() < freqder[chrom[ic].pop]) ? 1 : 0  ;
	config[newallele][chrom[ic].pop] += 1 ;
	if( (is < selspot) && ( (is+len) >= selspot ) ){
	  chrom[nchrom-2].allele = chrom[ic].allele ;
	  chrom[ic].allele = newallele ;
	  chrom[nchrom-1].allele = newallele ;
	}
	else {
	  chrom[nchrom-2].allele = newallele ;
	}
	ca( nsam,nsites, ic,  nchrom-1);
	return(ic);	

}

	int
xover(int nsam,int ic, int is)
{
	struct seg *pseg, *pseg2;
	int i,  lsg, lsgm1, newsg,  jseg, k,  in ;
	double ran1(), len ;


	pseg = chrom[ic].pseg ;
	lsg = chrom[ic].nseg ;
	len = (pseg + lsg -1)->end - pseg->beg ;
	cleft -= 1 - pow(pc,len) ;
   /* get seg # (jseg)  */

	for( jseg=0; is >= (pseg+jseg)->end ; jseg++) ;
	if( is >= (pseg+jseg)->beg ) in=1;
	else in=0;
	newsg = lsg - jseg ;

   /* copy last part of chrom to nchrom  */

	nchrom++;
	if( nchrom >= maxchr ) {
	    maxchr += 20 ;
	    chrom = (struct chromo *)realloc( chrom, (unsigned)(maxchr*sizeof(struct chromo))) ;
	    if( chrom == NULL ) perror( "malloc error. segtre2");
	    }
	if( !( pseg2 = chrom[nchrom-1].pseg = (struct seg *)calloc((unsigned)newsg,sizeof(struct seg)) ) )
		ERROR(" alloc error. re1");
	chrom[nchrom-1].nseg = newsg;
	chrom[nchrom-1].pop = chrom[ic].pop ;
	pseg2->end = (pseg+jseg)->end ;
	if( in ) {
		pseg2->beg = is + 1 ;
		(pseg+jseg)->end = is;
		}
	else pseg2->beg = (pseg+jseg)->beg ;
	pseg2->desc = (pseg+jseg)->desc ;
	for( k=1; k < newsg; k++ ) {
		(pseg2+k)->beg = (pseg+jseg+k)->beg;
		(pseg2+k)->end = (pseg+jseg+k)->end;
		(pseg2+k)->desc = (pseg+jseg+k)->desc;
		}

	lsg = chrom[ic].nseg = lsg-newsg + in ;
	lsgm1 = lsg - 1 ;
	nlinks -= pseg2->beg - (pseg+lsgm1)->end ;
	len = (pseg+lsgm1)->end - (pseg->beg) ;
	cleft += 1.0 - pow( pc, len) ;
	len = (pseg2 + newsg-1)->end - pseg2->beg ;
	cleft += 1.0 - pow(pc, len) ;
if( !(chrom[ic].pseg = 
     (struct seg *)realloc(chrom[ic].pseg,(unsigned)(lsg*sizeof(struct seg)) )) )
		perror( " realloc error. re2");
	if( in ) {
		begs = pseg2->beg;
		for( i=0,k=0; (k<nsegs-1)&&(begs > seglst[seglst[i].next].beg-1);
		   i=seglst[i].next, k++) ;
		if( begs != seglst[i].beg ) {
						/* new tree  */

	   	   if( nsegs >= seglimit ) {  
	   	   	  seglimit += SEGINC ;
	   	      nnodes = (int *)realloc( nnodes,(unsigned)(sizeof(int)*seglimit)) ; 
	   	      if( nnodes == NULL) perror("realloc error. 1. segtre_mig.c");
	   	      seglst =
	   	      	 (struct segl *)realloc( seglst,(unsigned)(sizeof(struct segl)*seglimit));
	   	      if(seglst == NULL ) perror("realloc error. 2. segtre_mig.c");
	   	      /*  printf("seglimit: %d\n",seglimit);  */
	   	      } 
	   	   seglst[nsegs].next = seglst[i].next;
	   	   seglst[i].next = nsegs;
	   	   seglst[nsegs].beg = begs ;
		   if( !(seglst[nsegs].ptree = (struct node *)calloc((unsigned)(2*nsam), sizeof(struct
			 node)) )) perror("calloc error. re3.");
		   nnodes[nsegs] = nnodes[i];
		   ptree1 = seglst[i].ptree;
		   ptree2 = seglst[nsegs].ptree;
		   nsegs++ ;
		   for( k=0; k<=nnodes[i]; k++) {
		      (ptree2+k)->abv = (ptree1+k)->abv ;
		      (ptree2+k)->time = (ptree1+k)->time;
		      }
		   }
	}
	return(ic) ;
}

/***** common ancestor subroutine **********************
   Pick two chromosomes and merge them. Update trees if necessary. **/

	int
ca(nsam, nsites,c1,c2)
	int nsam,c1,c2;
	int  nsites;
{
	int yes1, yes2, seg1, seg2, seg ;
	int tseg, start, end, desc, k;
	struct seg *pseg;
	struct node *ptree;

	seg1=0;
	seg2=0;

	if( !(pseg = (struct seg *)calloc((unsigned)nsegs,sizeof(struct seg) ))) 
		perror("alloc error.ca1");

	tseg = -1 ;

	for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
		start = seglst[seg].beg;
		yes1 = isseg(start, c1, &seg1);
		yes2 = isseg(start, c2, &seg2);
		if( yes1 || yes2 ) {
			tseg++;
			(pseg+tseg)->beg=seglst[seg].beg;
			end = ( k< nsegs-1 ? seglst[seglst[seg].next].beg-1 : nsites-1 ) ;
			(pseg+tseg)->end = end ;

			if( yes1 && yes2 ) {
				nnodes[seg]++;
				if( nnodes[seg] >= (2*nsam-2) ) tseg--;
				else
					(pseg+tseg)->desc = nnodes[seg];
				ptree=seglst[seg].ptree;
				desc = (chrom[c1].pseg + seg1) ->desc;
				(ptree+desc)->abv = nnodes[seg];
				desc = (chrom[c2].pseg + seg2) -> desc;
				(ptree+desc)->abv = nnodes[seg];
				(ptree+nnodes[seg])->time = t;

				}
			else {
				(pseg+tseg)->desc = ( yes1 ?
				   (chrom[c1].pseg + seg1)->desc :
				  (chrom[c2].pseg + seg2)->desc);
				}
			}
		}
	nlinks -= links(c1);
	cleft -= 1.0 - pow(pc, (double)links(c1));
	free(chrom[c1].pseg) ;
	if( tseg < 0 ) {
		free(pseg) ;
		chrom[c1].pseg = chrom[nchrom-1].pseg;
		chrom[c1].nseg = chrom[nchrom-1].nseg;
		chrom[c1].pop = chrom[nchrom-1].pop ;
		chrom[c1].allele = chrom[nchrom-1].allele ;
		if( c2 == nchrom-1 ) c2 = c1;
		nchrom--;
		}
	else {
		if( !(pseg = (struct seg *)realloc(pseg,(unsigned)((tseg+1)*sizeof(struct seg)))))
			perror(" realloc error. ca1");
		chrom[c1].pseg = pseg;
		chrom[c1].nseg = tseg + 1 ;
		nlinks += links(c1);
	   	cleft += 1.0 - pow(pc, (double)links(c1));
		}
	nlinks -= links(c2);
	cleft -= 1.0 - pow(pc, (double)links(c2));
	free(chrom[c2].pseg) ;
	chrom[c2].pseg = chrom[nchrom-1].pseg;
	chrom[c2].nseg = chrom[nchrom-1].nseg;
	chrom[c2].pop = chrom[nchrom-1].pop ;
	chrom[c2].allele = chrom[nchrom-1].allele ;
	nchrom--;
	if(tseg<0) return( 2 );  /* decrease of nchrom is two */
	else return( 1 ) ;
}

/*** Isseg: Does chromosome c contain the segment on seglst which starts at
	    start? *psg is the segment of chrom[c] at which one is to begin 
	    looking.  **/

	int
isseg(start, c, psg)
	int start, c, *psg;
{
	int ns;
	struct seg *pseg;

	ns = chrom[c].nseg;
	pseg = chrom[c].pseg;

/*  changed order of test conditions in following line on 6 Dec 2004 */
	for(  ; ((*psg) < ns ) && ( (pseg+(*psg))->beg <= start ) ; ++(*psg) )
		if( (pseg+(*psg))->end >= start ) return(1);
	return(0);
}
	


	void
	pick2_chrom(pop,config,allele,pc1,pc2)
	     int pop, *pc1, *pc2, config[],allele;
{
	int c1, c2, cs,cb,i, count;
	int pick2(int, int *, int *);
	pick2(config[pop],&c1,&c2);
	cs = (c1>c2) ? c2 : c1;
	cb = (c1>c2) ? c1 : c2 ;
		
	i=count=0;
	for(;;){
	  while( (chrom[i].pop != pop ) || ( chrom[i].allele != allele) ) i++;
		if( count == cs ) break;
		count++;
		i++;
		}
	*pc1 = i;
	i++;
	count++;
	for(;;){
	  while( (chrom[i].pop != pop) || (chrom[i].allele != allele)  ) i++;
		if( count == cb ) break;
		count++;
		i++;
		}
	*pc2 = i ;
}
	
	

/****  links(c): returns the number of links between beginning and end of chrom **/

	int
links(c)
	int c;
{
	int ns;

	ns = chrom[c].nseg - 1 ;

	return( (chrom[c].pseg + ns)->end - (chrom[c].pseg)->beg);
}

