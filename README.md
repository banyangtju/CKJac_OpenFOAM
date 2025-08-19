# CKJac installation and usage:
   The LoadBalancedChemistryModel is derived from DLBFoam (https://github.com/Aalto-CFD/DLBFoam), and integrates the recently developed **CKJac** (for more details in https://www.sciencedirect.com/science/article/abs/pii/S0010218024004978) chemistry model for computational speedup.

   The chemical model program, **CKJac**, provides the pre-computed and pre-saved fully analytic molar concentration-based and more spare Jacobian and chemical reaction sources, accelerating chemical reaction computation, which occupies 80% CPU time during CFD simulations. The speedup performance depends on chemical mechanisms, and more details can be found in our article.

   When using the CKJac chemistry model, a library called **libchemkin.so** should be generated in constant/chemkin. The libchemkin.so is produced based on chemkin-type reaction mechanisms （chem.bin). The Linux executable program JAC can be found in tutorials/chemFoam/CH4Air_GRIMech3.0/constant/chemkin. When using JAC, chem.bin should be provided in the same file directory and you can use the tutorials/chemFoam/CH4Air_GRIMech3.0/constant/chemkin/ckinterp/ to produce chem.bin. Thus, the operating steps are:

```
cd tutorials/chemFoam/CH4Air_GRIMech3.0/constant/chemkin/ckinterp/
make clean;make
# replace files (chem.inp and therm.dat) with your mechanism
./ckinterp
cp -r ./chem.bin ../
./JAC
```
   Finally, the library libchemkin.so can be produced.
   
   The advantages of CKJac can be found in Fig. 1 and Fig. 2. After installing the library **libchemistryModel_LB_LB_pyJac_LB_CKJac.so** and **libchemkin.so**, you can use it in constant/chemistryProperties,

```
chemistryType
{
    solver          ode;
    method          LB_CKJac;//LB;//LB_pyJac;//
}
```

and in system/controlDict
```
libs
(
	"$FOAM_CASE/constant/chemkin/libchemkin.so"
	"libchemistryModel_LB_LB_pyJac_LB_CKJac.so"
//	"libODE_seulex_LAPACK_seulex_KLU.so"
//	"$FOAM_CASE/constant/foam_mech/libc_pyjac.so"
)
```

# Attention:
   1. The ODE solvers including seulex_KLU (https://www.sciencedirect.com/science/article/pii/S001021801630267X) for CKJac (https://www.sciencedirect.com/science/article/abs/pii/S0010218024004978) and seulex_LAPACK (https://github.com/Aalto-CFD/DLBFoam) for pyJac (https://www.sciencedirect.com/science/article/pii/S0010465517300462) are uploaded. 
   2. CKJac can identify various common reactions (such as pressure-dependent reactions with Lindemann, SRI, and Troe forms, and Landau-Teller Formulation reactions).
   3. When you are using this solver to publish paper, please kindly consider to cite following papers:

   [1] Yangyang Ban, Fan Zhang, Naiyuan Zhang, Shenghui Zhong, Jiajian Zhu, Yiqiang Pei, The improved performance of plasma assisted combustion (PAC) simulations using the fully analytical Jacobian, Combustion and Flame, 2024, 270, 113788.


:blush:Contributor1: Yangyang Ban, Tianjin University, banyang@tju.edu.cn; 


:blush:Contributor2: Fan Zhang, associate professor, Tianjin University, fanzhang_lund@tju.edu.cn; 


:blush:Contributor3: Shenghui Zhong, Tianmushan Laboratory, https://github.com/ZSHtju. 


<div align=center>
	
**Fig. 1 The sparsity of CKJac.**
	
![image](https://github.com/user-attachments/assets/3f202ddb-3050-4888-8d06-fc7e96d06ec2)


**Fig. 2 The computational CPU time comparison.**
<div align=center>

![image](https://github.com/user-attachments/assets/0d170cf7-4d10-4ded-9267-55c165c5cbcf)


