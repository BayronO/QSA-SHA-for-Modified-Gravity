# Quasi-Static and Sub-Horizon Approximations in Modified Gravity

Files used to do all the computations in . You will find the following files:

1) Horndeski_Theory.nb: a Mathematica file containing the linear cosmological perturbations for the Horndeski model. You can find the perturbation equations for $f(R)$ theories in the section "Examples". You will need the xPand package in the xAct suite.

2) QSA_SHA_f(R).nb: a Mathematica file where the QSA and the SHA are applied to $f(R)$ theories. You will find two sections: $i)$ the standard procedure is applied, and $ii)$ the new parametrization is applied. You will need the xPand package in the xAct suite. 

3) fDES.nb: a Mathematica file where all the computations related to the $f$DES model are performed. In order to do the plots, you will need the data-files with the namefile "fDES_*.txt" and the TexAct package in the xAct suite.

4) Hu_Sawicky.nb: a Mathematica file where all the computations related to the Hu-Sawicky model are performed. In order to do the plots, you will need the data-files with the namefile "HS_*.txt" and the TexAct package in the xAct suite.

5) Full Solver: In this folder you will find the full solver of the perturbation equations for the $f$DES and Hu-Sawicki models. The solver is based on a Fortran routine known as dverk. To run the code using a terminal window:
    i) Go to the Full Solver folder using cd command
   ii) make clean
  iii) make HS - or - make designer
   iv) ./HS - or - ./designer

The data-files "base_ 2018_plikHM_TTTEEE_lowl_lowE_lensing00_pk.dat" and "base_2018_plikHM_TTTEEE_lowl_lowE_lensing00_cl.dat" contain data for the matter power spectrum and CMB angular power spectrum in the ΛCDM model computed using CLASS. All the data-files needed can be found in the folder "Data", where you also can find the data for the $P(k)$ and the $TT$ CMB power spectrum of the $f$DES model computed using EFCLASS.

