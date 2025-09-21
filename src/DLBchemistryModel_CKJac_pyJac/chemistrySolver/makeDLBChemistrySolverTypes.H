/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/


#ifndef makeDLBChemistrySolverTypes_H
#define makeDLBChemistrySolverTypes_H

#include "chemistrySolver.H"

#include "LoadBalancedChemistryModel.H"
#include "LoadBalanced_pyJacChemistryModel.H"
#include "LoadBalanced_CKJacChemistryModel.H"

#include "noChemistrySolver.H"
#include "EulerImplicit.H"
#include "ode.H"
#include "ode_pyJac.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDLBChemistrySolverType(SS, Comp, Thermo)                              \
                                                                               \
    typedef SS<LoadBalancedChemistryModel<Comp, Thermo>> LB##SS##Comp##Thermo;         \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        LB##SS##Comp##Thermo,                                                      \
        (#SS"<" + word(LoadBalancedChemistryModel<Comp, Thermo>::typeName_()) + "<"\
        + word(Comp::typeName_()) + "," + Thermo::typeName() + ">>").c_str(),  \
        0                                                                      \
    );                                                                         \
                                                                               \
    BasicChemistryModel<Comp>::                                                \
        add##thermo##ConstructorToTable<LB##SS##Comp##Thermo>                      \
        add##LB##SS##Comp##Thermo##thermo##ConstructorTo##BasicChemistryModel##Comp\
##Table_; \
                                                                               \
    typedef SS<LoadBalanced_pyJacChemistryModel<Comp, Thermo>> LBpyJac##SS##Comp##Thermo;       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        LBpyJac##SS##Comp##Thermo,                                                \
        (#SS"<" + word(LoadBalanced_pyJacChemistryModel<Comp, Thermo>::typeName_()) + "<"    \
        + word(Comp::typeName_()) + "," + Thermo::typeName() + ">>").c_str(),  \
        0                                                                      \
    );                                                                         \
                                                                               \
    BasicChemistryModel<Comp>::                                                \
        add##thermo##ConstructorToTable<LBpyJac##SS##Comp##Thermo>                \
        add##LBpyJac##SS##Comp##Thermo##thermo##ConstructorTo##BasicChemistryModel\
##Comp##Table_; \
                                                                               \
    typedef SS<LoadBalanced_CKJacChemistryModel<Comp, Thermo>> LBCKJac##SS##Comp##Thermo;       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        LBCKJac##SS##Comp##Thermo,                                                \
        (#SS"<" + word(LoadBalanced_CKJacChemistryModel<Comp, Thermo>::typeName_()) + "<"    \
        + word(Comp::typeName_()) + "," + Thermo::typeName() + ">>").c_str(),  \
        0                                                                      \
    );                                                                         \
                                                                               \
    BasicChemistryModel<Comp>::                                                \
        add##thermo##ConstructorToTable<LBCKJac##SS##Comp##Thermo>                \
        add##LBCKJac##SS##Comp##Thermo##thermo##ConstructorTo##BasicChemistryModel\
##Comp##Table_; 

#define makeDLBChemistrySolverTypes(Comp, Thermo)                                 \
                                                                               \
    makeDLBChemistrySolverType                                                    \
    (                                                                          \
        noChemistrySolver,                                                     \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeDLBChemistrySolverType                                                    \
    (                                                                          \
        EulerImplicit,                                                         \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeDLBChemistrySolverType                                                    \
    (                                                                          \
        ode,                                                                   \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
																			   \
    makeDLBChemistrySolverType                                                    \
    (                                                                          \
        ode_pyJac,                                                                   \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
