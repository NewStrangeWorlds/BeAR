/************************************************************************
 * $Id: cdisort.h 2887 2013-03-11 09:19:28Z robert.buras $
 ************************************************************************/

/*
 *   Copyright (c) 2011 by Timothy E. Dowling
 *   
 *   This file is part of cdisort.
 *
 *   cdisort is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   cdisort is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with cdisort.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#ifndef __cdisort_macros_h
#define __cdisort_macros_h


#if defined (__cplusplus)
extern "C" { 
#endif


/*
 * Array shift macros
 * Using unit-offset shift macros to match Fortran version
 *
 * NOTE: ARRAY(iq,jq) is defined locally instead of here, because its size is different
 *       in different subroutines.
 */
#define A(i,j)           a[i-1+(j-1)*lda]
#define AA(j,k)          aa[j-1+(k-1)*ia]
#define ABD(i,j)         abd[i-1+(j-1)*lda]
#define ALBMED(iu)       out->albmed[iu-1]
#define AMB(iq,jq)       ab[iq-1+(jq-1)*(ds->nstr/2)].zero
#define APB(iq,jq)       ab[iq-1+(jq-1)*(ds->nstr/2)].one

#define B(iq)            b[iq-1]
#define BDR(iq,jq)       bdr[iq-1+(jq)*(ds->nstr/2)]
#define BEM(iq)          bem[iq-1]

#define CBAND(irow,ncol) cband[irow-1+(ncol-1)*(9*(ds->nstr/2)-2)]
#define CC(iq,jq)        cc[iq-1+(jq-1)*ds->nstr]
#define CH(lc)           ch[lc-1]
#define CHTAU(ls)        chtau[ls]
#define CMU(iq)          cmu[iq-1]
#define CWT(iq)          cwt[iq-1]

#define DFDT(lu)         out->rad[lu-1].dfdt
#define DIAG(i)          diag[i-1].on
#define DTAUC(lc)        ds->dtauc[lc-1]
#define DTAU_C(lc)       dtau_c[lc-1]
#define DTAUCPR(lc)      dtaucpr[lc-1]

#define EMU(iu)          emu[iu-1]
#define EVAL(j)          eval[j-1]
#define EVEC(j,k)        evec[j-1+(k-1)*ievec]
#define EVECC(iq,jq)     evecc[iq-1+(jq-1)*ds->nstr]
#define EXPBEA(lc)       expbea[lc]

#define FLDIR(lu)        fl[lu-1].zero
#define FLDN(lu)         fl[lu-1].one
#define FLUP(lu)         out->rad[lu-1].flup
#define FLYR(lc)         flyr[lc-1]

#define GC(iq,jq,lc)     gc[iq-1+(jq-1+(lc-1)*ds->nstr)*ds->nstr]
#define GENSRC(maz,lc,iq)  ds->gensrc[iq-1+(lc-1+maz*ds->nlyr)*ds->nstr]
#define GENSRCU(maz,lc,iu) ds->gensrcu[iu-1+(lc-1+maz*ds->nlyr)*ds->numu]
#define GG(lc)           gg[lc-1]
#define GGPRIM(lc)       ggprim[lc-1]
#define GL(k,lc)         gl[k+(lc-1)*(ds->nstr+1)]
#define GMU(k)           gmu[k-1]
#define GU(iu,iq,lc)     gu[iu-1+(iq-1+(lc-1)*ds->nstr)*ds->numu]
#define GWT(k)           gwt[k-1]

#define IERROR(i)        ierror[i-1]
#define IPVT(k)          ipvt[k-1]

#define KK(iq,lc)        kk[iq-1+(lc-1)*ds->nstr]

#define LAYRU(lu)        layru[lu-1]
#define LL(iq,lc)        ll[iq-1+(lc-1)*ds->nstr]

#define MU(i)            mu[i-1]

#define OMEGA(lyr)       omega[lyr-1]
#define OPRIM(lc)        oprim[lc-1]

#define PKAG(lc)         pkag[lc]
#define PKAGC(lc)        pkagc[lc-1]
#define PHASA(lc)        phasa[lc-1]
#define PHASE(lc)        phase[lc-1]
#define PHASM(lc)        phasm[lc-1]
#define PHAST(lc)        phast[lc-1]
#define PHI(j)           ds->phi[j-1]
#define PHIRAD(jp)       phirad[jp-1]
#define PMOM(k,lc)       ds->pmom[k+(lc-1)*(ds->nmom_nstr+1)]
#define PRNTU0(i)        prntu0[i-1]
#define PSI0(iq)         psi[iq-1].zero
#define PSI1(iq)         psi[iq-1].one

#define RFLDIR(lu)       out->rad[lu-1].rfldir
#define RFLDN(lu)        out->rad[lu-1].rfldn
#define RMU(iu,iq)       rmu[iu-1+(iq)*ds->numu]
#define RR(lc)           rr[lc-1]

#define SSALB(lc)        ds->ssalb[lc-1]
#define SUBD(i)          diag[i-1].sub
#define SUPERD(i)        diag[i-1].super
#define SX(i)            sx[i-1]
#define SY(i)            sy[i-1]

#define TAU(lc)          tau[lc]
#define TAUC(lc)         tauc[lc]
#define TAUCPR(lc)       taucpr[lc]
#define TEMPER(lc)       ds->temper[lc]
#define TRNMED(iu)       out->trnmed[iu-1]

#define U0C(iq,lu)       u0c[iq-1+(lu-1)*ds->nstr]
#define U0U(iu,lu)       out->u0u[iu-1+(lu-1)*ds->numu]
#define UAVG(lu)         out->rad[lu-1].uavg
#define UAVGDN(lu)       out->rad[lu-1].uavgdn
#define UAVGUP(lu)       out->rad[lu-1].uavgup
#define UAVGSO(lu)       out->rad[lu-1].uavgso
#define UMU(iu)          ds->umu[iu-1]
#define UTAU(lu)         ds->utau[lu-1]
#define UTAUPR(lu)       utaupr[lu-1]
#define UUM(iu,lu)       uum[iu-1+(lu-1)*ds->numu]
#define UU(iu,lu,j)      out->uu[iu-1+(lu-1+(j-1)*ds->ntau)*ds->numu]
#define OUT_UUM(iu,lu,j) out->uum[iu-1+(lu-1+(j)*ds->ntau)*ds->numu] /* No -i behind j as mazim starts at 0, aky */

#define WK(iq)           wk[iq-1]

#define XBA(lc)          xba[lc]
#define XB0(iq,lc)       xb[iq-1+(lc-1)*ds->nstr].zero
#define XB1(iq,lc)       xb[iq-1+(lc-1)*ds->nstr].one
#define XB_0D(lc)        ts[lc-1].xb_0d
#define XB_0U(lc)        ts[lc-1].xb_0u
#define XB_1D(lc)        ts[lc-1].xb_1d
#define XB_1U(lc)        ts[lc-1].xb_1u
#define XP_0(lc)         ts[lc-1].xp_0
#define XP_1(lc)         ts[lc-1].xp_1
#define XR0(lc)          xr[lc-1].zero
#define XR1(lc)          xr[lc-1].one

#define YB_0D(lc)        ts[lc-1].yb_0d
#define YB_0U(lc)        ts[lc-1].yb_0u
#define YB_1D(lc)        ts[lc-1].yb_1d
#define YB_1U(lc)        ts[lc-1].yb_1u
#define YLM(l,i)         ylm[l+(i-1)*(maxmu+1)]
#define YLM0(iq)         ylm0[iq]
#define YLMC(l,iq)       ylmc[l+(iq-1)*(ds->nstr+1)]
#define YLMU(l,iu)       ylmu[l+(iu-1)*(ds->nstr+1)]
#define YP_0D(lc)        ts[lc-1].yp_0d
#define YP_0U(lc)        ts[lc-1].yp_0u
#define YP_1D(lc)        ts[lc-1].yp_1d
#define YP_1U(lc)        ts[lc-1].yp_1u

#define Z(j)             z[j-1]
#define Z0(iu)           zee[iu-1].zero
#define Z1(iq)           zee[iq-1].one
#define Z0U(iu,lc)       zu[iu-1+(lc-1)*ds->numu].zero
#define Z1U(iu,lc)       zu[iu-1+(lc-1)*ds->numu].one
#define ZB0U(iu,lc)      zbu[iu-1+(lc-1)*ds->numu].zero
#define ZB1U(iu,lc)      zbu[iu-1+(lc-1)*ds->numu].one
#define ZBAU(iu,lc)      zbu[iu-1+(lc-1)*ds->numu].alpha
#define ZB_A(lc)         ts[lc-1].zb_a
#define ZBEAM(iu,lc)     zbeam[iu-1+(lc-1)*ds->numu]
#define ZBEAMA(lc)       zbeama[lc-1]
#define ZBEAM0(iq,lc)    zbeamsp[iq-1+(lc-1)*ds->nstr].zero
#define ZBEAM1(iq,lc)    zbeamsp[iq-1+(lc-1)*ds->nstr].one
#define ZBS0(iq)         zbs[iq-1].zero
#define ZBS1(iq)         zbs[iq-1].one
#define ZD(j)            zd[j]
#define ZJ(j)            zj[j-1]
#define ZJG(j)           zjg[j-1]
#define ZJU(j)           zju[j-1]
#define ZGU(iu,lc)       zgu[iu-1+(lc-1)*ds->numu]
#define ZP_A(lc)         ts[lc-1].zp_a
#define ZPLK0(iq,lc)     plk[iq-1+(lc-1)*ds->nstr].zero
#define ZPLK1(iq,lc)     plk[iq-1+(lc-1)*ds->nstr].one
#define ZZ(iq,lc)        zz[iq-1+(lc-1)*ds->nstr]
#define ZZG(iq,lc)       zzg[iq-1+(lc-1)*ds->nstr]

/* BDE stuff */
#define MUP(it)          mu_phase[it-1]
#define PHASR(lc)        phasr[lc-1]
#define PHAS2(it,lc)     phas2[it-1+(lc-1)*nphase]
#define DSPHASE(it,lc)   ds->phase[it-1+(lc-1)*ds->nphase]
#define F_PHAS2_ABS(it)  f_phas2_abs[it-1]
#define MU_EQ(i,lu)      mu_eq[i-1+(lu-1)*nf]
#define NEG_PHAS(i,lu)   neg_phas[i-1+(lu-1)*nf]
#define NORM_PHAS(lu)    norm_phas[lu-1]

/* setout.f, inter.f stuff */
#define SDTAUC(i)        sdtauc[i-1]
#define SUTAU(i)         sutau[i-1]
#define ZOUT(i)          zout[i-1]
#define TAUINT(i)        tauint[i-1]
#define XARR(i)          xarr[i-1]
#define YARR(i)          yarr[i-1]

/*
 * Logical
 */
#define TRUE  1
#define FALSE 0

#define FIRST_IPHAS          1
#define ISOTROPIC            1
#define RAYLEIGH             2
#define HENYEY_GREENSTEIN    3
#define HAZE_GARCIA_SIEWERT  4
#define CLOUD_GARCIA_SIEWERT 5
#define LAST_IPHAS           5

#define GENERAL_BC 0
#define SPECIAL_BC 1

#define TOP_ILLUM 1
#define BOT_ILLUM 2

#define DS_WARNING 0
#define DS_ERROR   1

#define VERBOSE 0
#define QUIET   1

#define BRDF_NONE   0
#define BRDF_RPV    1  /* don't change these numbers as they are */
#define BRDF_CAM    2  /* used by Fortran code which of course   */
#define BRDF_AMB    3  /* has no access to this header file      */
#define BRDF_HAPKE  4

/*defined for new option names brdf_cam for cox_and_munk_sal,pcl,u10,uphi*/
#define BRDF_CAM_NN   4
#define BRDF_CAM_SAL  0 
#define BRDF_CAM_PCL  1
#define BRDF_CAM_U10  2
#define BRDF_CAM_UPHI 3

/*
 * NMUG : Number of angle cosine quadrature points on (-1,1) for integrating bidirectional reflectivity
 *        to get directional emissivity (it is necessary to use a quadrature set distinct from the 
 *        computational angles, because the computational angles may not be dense enough---ds->nstr
 *        may be too small---to give an accurate approximation for the integration).
 */
#define NMUG 50

/*
 * Mathematical
 */
#if !defined(M_E)
#  define M_E         2.7182818284590452354
#  define M_LOG2E     1.4426950408889634074
#  define M_LOG10E    0.43429448190325182765
#  define M_LN2       0.69314718055994530942
#  define M_LN10      2.30258509299404568402
#  define M_PI        3.14159265358979323846
#  define M_PI_2      1.57079632679489661923
#  define M_PI_4      0.78539816339744830962
#  define M_1_PI      0.31830988618379067154
#  define M_2_PI      0.63661977236758134308
#  define M_2_SQRTPI  1.12837916709551257390
#  define M_SQRT2     1.41421356237309504880
#  define M_SQRT1_2   0.70710678118654752440
#endif

#define DEG (M_PI/180.)

#define SQR(x) ({ \
          const double _x = (double)(x); \
          _x*_x; })

#define MIN(x,y) ({ \
         const double _x = (double)(x); \
         const double _y = (double)(y); \
         _x < _y ? _x : _y; })

#define MAX(x,y) ({ \
         const double _x = (double)(x); \
         const double _y = (double)(y); \
         _x > _y ? _x : _y; })

#define LIMIT_RANGE(min,x,max) ({ \
         const double _min = (double)(min); \
         const double _x   = (double)(x);   \
         const double _max = (double)(max); \
         _x < _min ? _min : ( _x > _max ? _max : _x ); })

#define IMIN(i,j) ({ \
         const int _i = (int)(i); \
         const int _j = (int)(j); \
         _i < _j ? _i : _j; })

#define IMAX(i,j) ({ \
         const int _i = (int)(i); \
         const int _j = (int)(j); \
         _i > _j ? _i : _j; })

#define F77_SIGN(a,b) ((b) >= 0. ? fabs(a) : -fabs(a))


/* * * * * * * * * * * * * * * * * * * * * * * end of cdisort.h * * * * * * * * * * * * * * * * * * */

#if defined (__cplusplus)
}
#endif


#endif

