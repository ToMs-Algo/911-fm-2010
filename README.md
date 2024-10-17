# Algorithm 911: Multiple-Precision Exponential Integral and Related Functions

> David M. Smith. 2011.  
> Algorithm 911: Multiple-Precision Exponential Integral and Related Functions.  
> ACM Trans. Math. Softw. 37, 4, Article 46 (February 2011), 16 pages.  
> https://doi.org/10.1145/1916461.1916470


## FM PACKAGE

See: [FM_User_Manual](Doc/FM_User_Manual.txt)

```for
!     FM 1.3                              David M. Smith                              6-21-2010
```

The routines in the file FM.f95 perform the low-level multiple precision arithmetic and functions
on three kinds of numbers.

FM routines handle floating-point real multiple precision numbers,
IM routines handle integer multiple precision numbers, and
ZM routines handle floating-point complex multiple precision numbers.

References to FM numbers below mean the low-level array form of the number used by the routines
in FM.f95, and not the derived type (fm) numbers handled by the FMZM module.  Logically, both
may refer to the same multiple precision number, but the syntax for dealing with the two types
of objects is different.  The same is true of references to IM numbers and ZM numbers below.

These are the basic routines for the FM package, and the expectation is that the user will not
call these routines directly.  The typical usage is for a program to declare multiple precision
variables with the three derived types defined in module FMZM in file FMZM90.f95.  Then that
module provides the interface between the user's program and the routines in this file.  See the
documentation in the sections above for advice on using the FMZM module.  The information below
is intended as a technical reference on the inner workings of FM, and most FM users should not
need to study it.


### LIST OF ROUTINES

#### Summary

```
--------------   Multiple precision versions of Fortran operations and functions   -----------------


   =
   +
   -
   *
   /
   **
   ==
   /=
   <
   <=
   >
   >=
   ABS          real    integer    complex
   ACOS         real               complex
   AIMAG                           complex
   AINT         real               complex
   ANINT        real               complex
   ASIN         real               complex
   ATAN         real               complex
   ATAN2        real
   BTEST                integer
   CEILING      real    integer    complex
   CMPLX        real    integer
   CONJG                           complex
   COS          real               complex
   COSH         real               complex
   DBLE         real    integer    complex
   DIGITS       real    integer    complex
   DIM          real    integer
   DINT         real               complex
   EPSILON      real
   EXP          real               complex
   EXPONENT     real
   FLOOR        real    integer    complex
   FRACTION     real               complex
   HUGE         real    integer    complex
   INT          real    integer    complex
   LOG          real               complex
   LOG10        real               complex
   MAX          real    integer
   MAXEXPONENT  real
   MIN          real    integer
   MINEXPONENT  real
   MOD          real    integer
   MODULO       real    integer
   NEAREST      real
   NINT         real    integer    complex
   PRECISION    real               complex
   RADIX        real    integer    complex
   RANGE        real    integer    complex
   REAL         real    integer    complex
   RRSPACING    real
   SCALE        real               complex
   SETEXPONENT  real
   SIGN         real    integer
   SIN          real               complex
   SINH         real               complex
   SPACING      real
   SQRT         real               complex
   TAN          real               complex
   TANH         real               complex
   TINY         real    integer    complex


-----------------------------   Conversion and inquiry functions   ---------------------------------


   TO_FM        real    integer    complex    string    other
   TO_IM        real    integer    complex    string    other
   TO_ZM        real    integer    complex    string    other
   TO_INT       real    integer    complex
   TO_SP        real    integer    complex
   TO_DP        real    integer    complex
   TO_SPZ       real    integer    complex
   TO_DPZ       real    integer    complex
   IS_OVERFLOW  real    integer    complex
   IS_UNDERFLOW real    integer    complex
   IS_UNKNOWN   real    integer    complex


-----------------------------------   Formatting functions   ---------------------------------------


   FM_FORMAT    real
   IM_FORMAT            integer
   ZM_FORMAT                       complex


------------------------------------   Integer functions   -----------------------------------------


   FACTORIAL            integer
   GCD                  integer
   MULTIPLY_MOD         integer
   POWER_MOD            integer


------------------------------------   Special functions   -----------------------------------------


   BERNOULLI(N)              real
   BESSEL_J(N,X)             real
   BESSEL_Y(N,X)             real
   BETA(A,B)                 real
   BINOMIAL(A,B)             real
   COS_INTEGRAL(X)           real
   COSH_INTEGRAL(X)          real
   ERF(X)                    real
   ERFC(X)                   real
   EXP_INTEGRAL_EI(X)        real
   EXP_INTEGRAL_EN(N,X)      real
   FACTORIAL(X)              real
   FRESNEL_C(X)              real
   FRESNEL_S(X)              real
   GAMMA(X)                  real
   INCOMPLETE_BETA(X,A,B)    real
   INCOMPLETE_GAMMA1(A,X)    real
   INCOMPLETE_GAMMA2(A,X)    real
   LOG_ERFC(X)               real
   LOG_GAMMA(X)              real
   LOG_INTEGRAL(X)           real
   POLYGAMMA(N,X)            real
   POCHHAMMER(X,N)           real
   PSI(X)                    real
   SIN_INTEGRAL(X)           real
   SINH_INTEGRAL(X)          real
```

Several of these functions are described in more detail below.

#### FM - multiple precision real numbers

> 9.  LIST OF ROUTINES

First are the routines that deal with multiple precision real numbers.  All of these are
subroutines except logical function FMCOMPARE.

MA, MB, MC refer to FM format numbers (i.e., integers as opposed to the derived types with
integer components that are defined in file FMZM90.f95)

In Fortran-90 and later versions of the Fortran standard, it is potentially unsafe to use the
same variable more than once in the calling sequence.  The operation MA = MA + MB should not be
written as

      CALL FMADD(MA,MB,MA)

since the code for the subroutine will not know that the first and third arguments are the same,
and some code optimizations under the assumption that all three arguments are different could
cause errors.

One solution is to use a third array and then put the result back in MA:

      CALL FMADD(MA,MB,MC)
      CALL FMEQ(MC,MA)

When the first call is doing one of the "fast" operations like addition, the extra call to move
the result back to MA can cause a noticeable loss in efficiency.  To avoid this, separate
routines are provided for the basic arithmetic operations when the result is to be returned in
the same array as one of the inputs.

A routine name with a suffix of  "_R1" returns the result in the first input array, and a suffix
of "_R2" returns the result in the second input array.  The example above would then be:

      CALL FMADD_R1(MA,MB)

These routines each have one less argument than the original version, since the output is
re-directed to one of the inputs.  The result array should not be the same as any input array
when the original version of the routine is used.

The routines that can be used this way are listed below.  For others, like

      CALL FMEXP(MA,MA)

the relative cost of doing an extra copy is small.  This one should become

      CALL FMEXP(MA,MB)
      CALL FMEQ(MB,MA)

When the derived-type interface is used, as in

      TYPE (FM), SAVE :: A, B
      ...
      A = A + B

there is no problem putting the result back into A, since the interface routine creates a
temporary scratch array for the result of A + B.

For each of these routines there is also a version available for which the argument list is
the same but all FM numbers are in packed format.  The routines using packed numbers have the
same names except 'FM' is replaced by 'FP' at the start of each name.

Some of the routine names were restricted to 6 characters in earlier versions of FM.  The old
names have been retained for compatibility, but new names that are longer and more readable
have been added.  For example, the old routine FMCSSN can now also be called as FMCOS_SIN.
Both old and new names are listed below.

```for
FMABS(MA,MB)              MB = ABS(MA)

FMACOS(MA,MB)             MB = ACOS(MA)

FMADD(MA,MB,MC)           MC = MA + MB

FMADD_R1(MA,MB)           MA = MA + MB

FMADD_R2(MA,MB)           MB = MA + MB

FMADDI(MA,IVAL)           MA = MA + IVAL   Increment an FM number by a one word integer.
                                           Note this call does not have an "MB" result
                                           like FMDIVI and FMMPYI.

FMASIN(MA,MB)             MB = ASIN(MA)

FMATAN(MA,MB)             MB = ATAN(MA)

FMATAN2(MA,MB,MC)         MC = ATAN2(MA,MB)     < old name: FMATN2 >

FMBIG(MA)                 MA = Biggest FM number less than overflow.

FMCHANGEBASE(MA,MB,NEW_MBASE,NEW_NDIG)
                          MB is returned with the base NEW_MBASE and precision NEW_NDIG
                             representation MA, where MA is given in the current base (MBASE)
                             and precision (NDIG).  This routine is primarily meant to be used
                             for input and output conversion when a base is being used that is
                             not a power of ten.

FMCOMPARE(MA,LREL,MB)     Logical comparison of MA and MB.     < old name: FMCOMP >
                          LREL is a character(2) value identifying which of the six comparisons
                               is to be made.
                          Example:  IF (FMCOMPARE(MA,'>=',MB)) ...
                          Also can be:  IF (FMCOMPARE(MA,'GE',MB)) ...
                          character(1) is ok:  IF (FMCOMPARE(MA,'>',MB)) ...

FMCONS                    Set several saved constants that depend on MBASE, the base being used.
                          FMCONS should be called immediately after changing MBASE.

FMCOS(MA,MB)              MB = COS(MA)

FMCOS_SIN(MA,MB,MC)       MB = COS(MA),  MC = SIN(MA).     < old name: FMCSSN >
                               Faster than making two separate calls.

FMCOSH(MA,MB)             MB = COSH(MA)

FMCOSH_SINH(MA,MB,MC)     MB = COSH(MA),  MC = SINH(MA).     < old name: FMCHSH >
                               Faster than making two separate calls.

FMDIG(NSTACK,KST)         Find a set of precisions to use during Newton iteration for finding a
                          simple root starting with about double precision accuracy.

FMDIM(MA,MB,MC)           MC = DIM(MA,MB)

FMDIV(MA,MB,MC)           MC = MA / MB

FMDIV_R1(MA,MB)           MA = MA / MB

FMDIV_R2(MA,MB)           MB = MA / MB

FMDIVI(MA,IVAL,MB)        MB = MA/IVAL   IVAL is a one word integer.

FMDIVI_R1(MA,IVAL)        MA = MA/IVAL

FMDP2M(X,MA)              MA = X    Convert from double precision to FM.

FMDPM(X,MA)               MA = X    Convert from double precision to FM.
                                    Faster than FMDP2M, but MA agrees with X only to D.P.
                                    accuracy.  See the comments in the two routines.

FMEQ(MA,MB)               MB = MA   Both have precision NDIG.
                                    This is the version to use for standard  B = A  statements.

FMEQU(MA,MB,NA,NB)        MB = MA   Version for changing precision.
                                    MA has NA digits (i.e., MA was computed using NDIG = NA), and
                                    MB will be defined having NB digits.
                                    MB is rounded if NB < NA
                                    MB is zero-padded if NB > NA

FMEXP(MA,MB)              MB = EXP(MA)

FMFLAG(K)                 K = KFLAG  get the value of the FM condition flag -- stored in the
                                     internal FM variable KFLAG in module FMVALS.

FMFORM(FORM,MA,STRING)    MA is converted to a character string using format FORM and returned in
                             STRING.  FORM can represent I, F, E, or ES formats.  Example:
                             CALL FMFORM('F60.40',MA,STRING)

FMFPRINT(FORM,MA)         Print MA on unit KW using FORM format.     < old name: FMFPRT >

FMI2M(IVAL,MA)            MA = IVAL   Convert from one word integer to FM.

FMINP(LINE,MA,LA,LB)      MA = LINE   Input conversion.
                                      Convert LINE(LA) through LINE(LB) from characters to FM.

FMINT(MA,MB)              MB = INT(MA)    Integer part of MA.

FMIPOWER(MA,IVAL,MB)      MB = MA**IVAL   Raise an FM number to a one word integer power.
                                          < old name: FMIPWR >

FMLOG10(MA,MB)            MB = LOG10(MA)     < old name: FMLG10 >

FMLN(MA,MB)               MB = LOG(MA)

FMLNI(IVAL,MA)            MA = LOG(IVAL)   Natural log of a one word integer.

FMM2DP(MA,X)              X  = MA     Convert from FM to double precision.

FMM2I(MA,IVAL)            IVAL = MA   Convert from FM to integer.

FMM2SP(MA,X)              X  = MA     Convert from FM to single precision.

FMMAX(MA,MB,MC)           MC = MAX(MA,MB)

FMMIN(MA,MB,MC)           MC = MIN(MA,MB)

FMMOD(MA,MB,MC)           MC = MA mod MB

FMMPY(MA,MB,MC)           MC = MA * MB

FMMPY_R1(MA,MB)           MA = MA * MB

FMMPY_R2(MA,MB)           MB = MA * MB

FMMPYI(MA,IVAL,MB)        MB = MA*IVAL    Multiply by a one word integer.

FMMPYI_R1(MA,IVAL)        MA = MA*IVAL

FMNINT(MA,MB)             MB = NINT(MA)   Nearest FM integer.

FMOUT(MA,LINE,LB)         LINE = MA   Convert from FM to character.
                                      LINE is a character array of length LB.

FMPI(MA)                  MA = pi

FMPRINT(MA)               Print MA on unit KW using current format.     < old name: FMPRNT >

FMPOWER(MA,MB,MC)         MC = MA**MB     < old name: FMPWR >

FM_RANDOM_NUMBER(X)       X is returned as a double precision random number, uniformly
                          distributed on the open interval (0,1).  It is a high-quality,
                          long-period generator based on 49-digit prime numbers.
                          Note that X is double precision, unlike the similar Fortran intrinsic
                          random number routine, which returns a single-precision result.
                          A default initial seed is used if FM_RANDOM_NUMBER is called without
                          calling FM_RANDOM_SEED_PUT first.  See the comments in section 11 below
                          and also those in the routine for more details.

FM_RANDOM_SEED_GET(SEED)  returns the seven integers SEED(1) through SEED(7) as the current seed
                          for the FM_RANDOM_NUMBER generator.

FM_RANDOM_SEED_PUT(SEED)  initializes the FM_RANDOM_NUMBER generator using the seven integers
                          SEED(1) through SEED(7). These get and put functions are slower than
                          FM_RANDOM_NUMBER, so FM_RANDOM_NUMBER should be called many times
                          between FM_RANDOM_SEED_PUT calls.  Also, some generators that used a
                          9-digit modulus have failed randomness tests when used with only a few
                          numbers being generated between calls to re-start with a new seed.

FM_RANDOM_SEED_SIZE(SIZE) returns integer SIZE as the size of the SEED array used by the
                          FM_RANDOM_NUMBER generator.  Currently, SIZE = 7.

FMRATIONAL_POWER(MA,K,J,MB)
                          MB = MA**(K/J)  Rational power.     < old name: FMRPWR >
                          Faster than FMPOWER for functions like the cube root.

FMREAD(KREAD,MA)          MA   is returned after reading one (possibly multi-line) FM number
                               on unit KREAD.  This routine reads numbers written by FMWRITE.

FMSET(NPREC)              Set the internal FM variables so that the precision is at least NPREC
                          base 10 digits plus three base 10 guard digits.

FMSETVAR(STRING)          Define a new value for one of the internal FM variables in module
                          FMVALS that controls one of the FM options.  STRING has the form
                                variable = value.
                          Example:  To change the screen width for FM output:
                                CALL FMSETVAR(' KSWIDE = 120 ')
                          The variables that can be changed and the options they control are
                          listed in sections 2 through 6 above.  Only one variable can be set
                          per call.  The variable name in STRING must have no embedded blanks.
                          The value part of STRING can be in any numerical format, except in
                          the case of variable CMCHAR, which is character type.  To set CMCHAR
                          to 'E', don't use any quotes in STRING:
                                CALL FMSETVAR(' CMCHAR = E ')

FMSIGN(MA,MB,MC)          MC = SIGN(MA,MB)   Returns the absolute value of MA times the sign
                                             of MB.

FMSIN(MA,MB)              MB = SIN(MA)

FMSINH(MA,MB)             MB = SINH(MA)

FMSP2M(X,MA)              MA = X   Convert from single precision to FM.

FMSQR(MA,MB)              MB = MA * MA   Faster than FMMPY.

FMSQR_R1(MA)              MA = MA * MA

FMSQRT(MA,MB)             MB = SQRT(MA)

FMSQRT_R1(MA)             MA = SQRT(MA)

FMST2M(STRING,MA)         MA = STRING
                               Convert from character string to FM.  STRING may be in any
                               numerical format.  FMST2M is often more convenient than FMINP,
                               which converts an array of character(1) values.  Example:
                                     CALL FMST2M('123.4',MA)

FMSUB(MA,MB,MC)           MC = MA - MB

FMSUB_R1(MA,MB)           MA = MA - MB

FMSUB_R2(MA,MB)           MB = MA - MB

FMTAN(MA,MB)              MB = TAN(MA)

FMTANH(MA,MB)             MB = TANH(MA)

FMTINY(MA)                MA = Smallest positive FM number greater than underflow.

FMULP(MA,MB)              MB = One Unit in the Last Place of MA.  For positive MA this is the
                               same as the Fortran function SPACING, but MB < 0 if MA < 0.
                               Examples:  If MBASE = 10 and NDIG = 30, then ulp(1.0) = 1.0E-29,
                                          ulp(-4.5E+67) = -1.0E+38.

FMVARS                    Write the current values of the internal FM variables on unit KW.

FMWRITE(KWRITE,MA)        Write MA on unit KWRITE.     < old name: FMWRIT >
                          Multi-line numbers will have '&' as the last nonblank character on all
                          but the last line.  These numbers can then be read easily using FMREAD.
```

#### FM - special functions

These are the available mathematical special functions.

```for
FMBERNOULLI(N,MA)         MA = B(N)      Nth Bernoulli number

FMBESJ(N,MA,MB)           MB = J(N,MA)   Bessel function of the first kind.

FMBESY(N,MA,MB)           MB = Y(N,MA)   Bessel function of the second kind.

FMBETA(MA,MB,MC)          MC = Beta(MA,MB)

FMC(MA,MB)                MB = C(MA)     Fresnel Cosine Integral

FMCHI(MA,MB)              MB = Chi(MA)   Hyperbolic Cosine Integral

FMCI(MA,MB)               MB = Ci(MA)    Cosine Integral

FMCOMB(MA,MB,MC)          MC = Combination MA choose MB  (Binomial coefficient)

FMEI(MA,MB)               MB = Ei(MA)    Exponential Integral

FMEN(N,MA,MB)             MB = E(N,MA)   Exponential Integral E_n

FMERF(MA,MB)              MB = Erf(MA)   Error function

FMERFC(MA,MB)             MB = Erfc(MA)  Complimentary Error function

FMEULER(MA)               MA = Euler's constant ( 0.5772156649... )     < old name: FMEULR >

FMFACT(MA,MB)             MB = MA Factorial  (Gamma(MA+1))

FMGAM(MA,MB)              MB = Gamma(MA)

FMIBTA(MX,MA,MB,MC)       MC = Incomplete Beta(MX,MA,MB)

FMIGM1(MA,MB,MC)          MC = Incomplete Gamma(MA,MB).  Lower case Gamma(a,x)

FMIGM2(MA,MB,MC)          MC = Incomplete Gamma(MA,MB).  Upper case Gamma(a,x)

FMLERC(MA,MB)             MB = Ln(Erfc(MA))  Log Erfc

FMLI(MA,MB)               MB = Li(MA)    Logarithmic Integral

FMLNGM(MA,MB)             MB = Ln(Gamma(MA))

FMPGAM(N,MA,MB)           MB = Polygamma(N,MA)  (Nth derivative of Psi)

FMPOCH(MA,N,MB)           MB = MA*(MA+1)*(MA+2)*...*(MA+N-1)  (Pochhammer)

FMPSI(MA,MB)              MB = Psi(MA)   (Derivative of Ln(Gamma(MA))

FMS(MA,MB)                MB = S(MA)     Fresnel Sine Integral

FMSHI(MA,MB)              MB = Shi(MA)   Hyperbolic Sine Integral

FMSI(MA,MB)               MB = Si(MA)    Sine Integral
```

#### IM - multiple precision integer numbers

These are the routines that deal with multiple precision integer numbers.
All are subroutines except logical function IMCOMPARE.  MA, MB, MC refer to IM format numbers.
In each case the version of the routine to handle packed IM numbers has the same name, with
'IM' replaced by 'IP'.

```for
IMABS(MA,MB)              MB = ABS(MA)

IMADD(MA,MB,MC)           MC = MA + MB

IMBIG(MA)                 MA = 10**(10**6).
                               Larger IM numbers can be obtained, but setting MA to the largest
                               possible value would leave no room for any other numbers.

IMCOMPARE(MA,LREL,MB)     Logical comparison of MA and MB.     < old name: IMCOMP >
                          LREL is a character(2) value identifying which of the six comparisons
                               is to be made.
                          Example:  IF (IMCOMPARE(MA,'GE',MB)) ...
                          Also can be:  IF (IMCOMPARE(MA,'>=',MB))
                          character(1) is ok:  IF (IMCOMPARE(MA,'>',MB)) ...

IMDIM(MA,MB,MC)           MC = DIM(MA,MB)

IMDIV(MA,MB,MC)           MC = int(MA/MB)
                               Use IMDIVR if the remainder is also needed.

IMDIVI(MA,IVAL,MB)        MB = int(MA/IVAL)
                               IVAL is a one word integer.  Use IMDVIR to get the remainder also.

IMDIVR(MA,MB,MC,MD)       MC = int(MA/MB),   MD = MA mod MB
                               When both the quotient and remainder are needed, this routine is
                               twice as fast as calling both IMDIV and IMMOD.

IMDVIR(MA,IVAL,MB,IREM)   MB = int(MA/IVAL),   IREM = MA mod IVAL
                          IVAL and IREM are one word integers.

IMEQ(MA,MB)               MB = MA

IMFM2I(MAFM,MB)           MB = MAFM  Convert from real (FM) format to integer (IM) format.

IMFORM(FORM,MA,STRING)    MA is converted to a character string using format FORM and
                             returned in STRING.  FORM can represent I, F, E, or ES formats.
                             Example: CALL IMFORM('I70',MA,STRING)

IMFPRINT(FORM,MA)         Print MA on unit KW using FORM format.     < old name: IMFPRT >

IMGCD(MA,MB,MC)           MC = greatest common divisor of MA and MB.

IMI2FM(MA,MBFM)           MBFM = MA  Convert from integer (IM) format to real (FM) format.

IMI2M(IVAL,MA)            MA = IVAL   Convert from one word integer to IM.

IMINP(LINE,MA,LA,LB)      MA = LINE   Input conversion.
                                      Convert LINE(LA) through LINE(LB) from characters to IM.

IMM2DP(MA,X)              X  = MA     Convert from IM to double precision.

IMM2I(MA,IVAL)            IVAL = MA   Convert from IM to one word integer.

IMM2SP(MA,X)              X  = MA     Convert from IM to single precision.

IMMAX(MA,MB,MC)           MC = MAX(MA,MB)

IMMIN(MA,MB,MC)           MC = MIN(MA,MB)

IMMOD(MA,MB,MC)           MC = MA mod MB

IMMPY(MA,MB,MC)           MC = MA*MB

IMMPYI(MA,IVAL,MB)        MB = MA*IVAL    Multiply by a one word integer.

IMMPY_MOD(MA,MB,MC,MD)    MD = MA*MB mod MC     < old name: IMMPYM >
                               Slightly faster than calling IMMPY and IMMOD separately.

IMOUT(MA,LINE,LB)         LINE = MA   Convert from IM to character.
                                      LINE is a character array of length LB.

IMPOWER(MA,MB,MC)         MC = MA**MB     < old name: IMPWR >

IMPOWER_MOD(MA,MB,MC,MD)  MD = MA**MB mod MC     < old name: IMPMOD >

IMPRINT(MA)               Print MA on unit KW.     < old name: IMPRNT >

IMREAD(KREAD,MA)          MA   is returned after reading one (possibly multi-line)
                               IM number on unit KREAD.
                               This routine reads numbers written by IMWRITE.

IMSIGN(MA,MB,MC)          MC = SIGN(MA,MB)   Returns the absolute value of MA times the
                                             sign of MB.

IMSQR(MA,MB)              MB = MA*MA   Faster than IMMPY.

IMST2M(STRING,MA)         MA = STRING
                               Convert from character string to IM.
                               IMST2M is often more convenient than IMINP, which converts
                               an array of character(1) values.  Example:
                                    CALL IMST2M('12345678901',MA)

IMSUB(MA,MB,MC)           MC = MA - MB

IMWRITE(KWRITE,MA)        Write MA on unit KWRITE.
                          Multi-line numbers will have '&' as the last nonblank character on all
                          but the last line.  These numbers can then be read easily using IMREAD.
```

#### ZM - multiple precision complex numbers

These are the routines that deal with multiple precision complex numbers.
All are subroutines, and in each case the version of the routine to handle packed ZM numbers has
the same name, with 'ZM' replaced by 'ZP'.

```for
MA, MB, MC refer to ZM format complex numbers.
MAFM, MBFM, MCFM refer to FM format real numbers.
INTEG is a Fortran INTEGER variable.
ZVAL is a Fortran COMPLEX variable.

ZMABS(MA,MBFM)            MBFM = ABS(MA)    Result is real.

ZMACOS(MA,MB)             MB = ACOS(MA)

ZMADD(MA,MB,MC)           MC = MA + MB

ZMADDI(MA,INTEG)          MA = MA + INTEG  Increment an ZM number by a one word integer.
                                           Note this call does not have an "MB" result
                                           like ZMDIVI and ZMMPYI.

ZMARG(MA,MBFM)            MBFM = Argument(MA)    Result is real.

ZMASIN(MA,MB)             MB = ASIN(MA)

ZMATAN(MA,MB)             MB = ATAN(MA)

ZMCOMPLEX(MAFM,MBFM,MC)   MC = CMPLX(MAFM,MBFM)     < old name: ZMCMPX >

ZMCONJUGATE(MA,MB)        MB = CONJG(MA)     < old name: ZMCONJ >

ZMCOS(MA,MB)              MB = COS(MA)

ZMCOS_SIN(MA,MB,MC)       MB = COS(MA),  MC = SIN(MA).     < old name: ZMCSSN >
                               Faster than 2 calls.

ZMCOSH(MA,MB)             MB = COSH(MA)

ZMCOSH_SINH(MA,MB,MC)     MB = COSH(MA),  MC = SINH(MA).     < old name: ZMCHSH >
                               Faster than 2 calls.

ZMDIV(MA,MB,MC)           MC = MA / MB

ZMDIVI(MA,INTEG,MB)       MB = MA / INTEG

ZMEQ(MA,MB)               MB = MA

ZMEQU(MA,MB,NDA,NDB)      MB = MA    Version for changing precision.
                                     (NDA and NDB are as in FMEQU)

ZMEXP(MA,MB)              MB = EXP(MA)

ZMFORM(FORM1,FORM2,MA,STRING)
                          STRING = MA
                          MA is converted to a character string using format FORM1 for the real
                          part and FORM2 for the imaginary part.  The result is returned in
                          STRING.  FORM1 and FORM2 can represent I, F, E, or ES formats.
                          Example:
                                CALL ZMFORM('F20.10','F15.10',MA,STRING)

ZMFPRINT(FORM1,FORM2,MA)  Print MA on unit KW using formats FORM1 and FORM2.
                          < old name: ZMFPRT >

ZMI2M(INTEG,MA)           MA = CMPLX(INTEG,0)

ZM2I2M(INTEG1,INTEG2,MA)  MA = CMPLX(INTEG1,INTEG2)

ZMIMAG(MA,MBFM)           MBFM = IMAG(MA)    Imaginary part.

ZMINP(LINE,MA,LA,LB)      MA = LINE   Input conversion.
                               Convert LINE(LA) through LINE(LB) from characters to ZM.
                               LINE is a character array of length at least LB.

ZMINT(MA,MB)              MB = INT(MA)       Integer part of both Real and Imaginary parts of MA.

ZMIPOWER(MA,INTEG,MB)     MB = MA ** INTEG   Integer power function.     < old name: ZMIPWR >

ZMLOG10(MA,MB)            MB = LOG10(MA)     < old name: ZMLG10 >

ZMLN(MA,MB)               MB = LOG(MA)

ZMM2I(MA,INTEG)           INTEG = INT(REAL(MA))

ZMM2Z(MA,ZVAL)            ZVAL = MA

ZMMPY(MA,MB,MC)           MC = MA * MB

ZMMPYI(MA,INTEG,MB)       MB = MA * INTEG

ZMNINT(MA,MB)             MB = NINT(MA)   Nearest integer of both Real and Imaginary.

ZMOUT(MA,LINE,LB,LAST1,LAST2)
                          LINE = MA
                          Convert from FM to character.
                          LINE  is the returned character(1) array.
                          LB    is the dimensioned size of LINE.
                          LAST1 is returned as the position in LINE of the last character
                                of REAL(MA)
                          LAST2 is returned as the position in LINE of the last character
                                of AIMAG(MA)

ZMPOWER(MA,MB,MC)         MC = MA ** MB     < old name: ZMPWR >

ZMPRINT(MA)               Print MA on unit KW using current format.     < old name: ZMPRNT >

ZMRATIONAL_POWER(MA,IVAL,JVAL,MB)
                          MB = MA ** (IVAL/JVAL)     < old name: ZMRPWR >

ZMREAD(KREAD,MA)          MA   is returned after reading one (possibly multi-line) ZM number on
                               unit KREAD.  This routine reads numbers written by ZMWRITE.

ZMREAL(MA,MBFM)           MBFM = REAL(MA)    Real part.

ZMSET(NPREC)              Set precision to the equivalent of a few more than NPREC base 10
                          digits.  This is now the same as FMSET, but is retained for
                          compatibility with earlier versions of the package.

ZMSIN(MA,MB)              MB = SIN(MA)

ZMSINH(MA,MB)             MB = SINH(MA)

ZMSQR(MA,MB)              MB = MA*MA    Faster than ZMMPY.

ZMSQRT(MA,MB)             MB = SQRT(MA)

ZMST2M(STRING,MA)         MA = STRING
                               Convert from character string to ZM.  ZMST2M is often more
                               convenient than ZMINP, which converts an array of character(1)
                               values.  Example:
                                     CALL ZMST2M('123.4+5.67i',MA).

ZMSUB(MA,MB,MC)           MC = MA - MB

ZMTAN(MA,MB)              MB = TAN(MA)

ZMTANH(MA,MB)             MB = TANH(MA)

ZMWRITE(KWRITE,MA)        Write MA on unit KWRITE.  Multi-line numbers are formatted for
                          automatic reading with ZMREAD.     < old name: ZMWRIT >

ZMZ2M(ZVAL,MA)            MA = ZVAL
```


### NEW FOR VERSION 1.3

> 11. NEW FOR VERSION 1.3

The routines for the exponential integral function and related mathematical special functions
are new in version 1.3.  These routines are:
FMBESJ, FMBESY, FMC, FMCHI, FMCI, FMEI, FMEN, FMERF, FMERFC, FMLERC, FMLI, FMS, FMSHI, FMSI.

Some of the routines were moved between files FM.f95 and FMZM90.f95 so that now all routines
using the module FMZM (in file FMZM90.f95) for multiple precision derived types and operator
overloading are located in FMZM90.f95.  This means that programs not using derived types can
skip compiling and/or linking FMZM90.f95.

The array function DOTPRODUCT in FMZM has been re-named DOT_PRODUCT to agree with the Fortran
standard.  For type ZM complex arguments its definition has been changed to agree with the
Fortran intrinsic function.  When X and Y are complex, DOT_PRODUCT(X,Y) is not just the sum of
the products of the corresponding array elements, as it is for types FM and IM.  For type ZM,
the formula is the sum of conjg(X(j)) * Y(j).  This definition is used so that the complex dot
product will be an inner product in the mathematical sense.

New routines have been added to module FMZM to provide array syntax for the three multiple
precision derived types.  This means statements like V = 1 and A = B + C now work when these
variables are vectors or matrices of multiple precision numbers.

One routine from FM 1.2 has been split into three routines in version 1.3.  The routine
FM_RANDOM_SEED from FM 1.2 has become three subroutines, so that the optional arguments and
the need for an explicit interface can be avoided.  See the three routines starting with
FM_RANDOM_SEED in the list above.  The same multiplicative congruential generator as before
is used, but the shuffling of those values has been removed, so that saving seeds and
re-starting the generator now works more like the standard Fortran random function.

Multiple precision variables were separate fixed-size arrays in previous versions.  Now they are
single integers that serve as index values to a single large array (MWK, defined in file
FMSAVE.f95) where the actual values are stored.  This often improves both efficiency and memory
utilization, since many compilers implemented the derived type operations using copy in and copy
out of the arguments for a given operation.  Copying entire arrays was slower, and there were
often memory leaks when the compiler automatically created temporary derived type objects while
evaluating derived type expressions.  The static arrays in previous versions also meant that
memory was wasted when only a few kinds of operations were used at high precision.  Now the
space needed by any unused operations never gets allocated.

Some new error checking is now done for the derived type multiple precision variables. Attempting
to use an undefined variable will cause an error message to be printed.

Much higher precision can be attained in version 1.3, since machines are faster and have more
memory.  To support higher precision, a routine for FFT-based multiplication has been included,
and when precision gets high enough, the algorithms for multiplication, division, squares, square
roots, etc., will switch to the FFT routine.

Binary splitting algorithms are used for the mathematical constants at high precision.  At the
time version 1.3 was released, computing a million digits of e, pi, or the logarithm of a small
integer took a few seconds, while a million digits of Euler's constant took a few minutes.

Perfect rounding is now done all the time.  In version 1.2 perfect rounding was an option, but
the default rounding could round the wrong direction once every few million operations, when the
exact result was very close to halfway between two adjacent representable numbers.
