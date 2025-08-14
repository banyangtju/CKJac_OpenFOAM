/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DLBFoam: Dynamic Load Balancing 
   \\    /   O peration     | for fast reactive simulations
    \\  /    A nd           | 
     \\/     M anipulation  | 2020, Aalto University, Finland
-------------------------------------------------------------------------------
License
    This file is part of DLBFoam library, derived from OpenFOAM.

    https://github.com/blttkgl/DLBFoam

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

Class
    Foam::LoadBalanced_CKJacChemistryModel

Description
    Extends loadBalancedChemistryModel model by adding analytical Jacobian.

SourceFiles
    LoadBalanced_CKJacChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef LoadBalanced_CKJacChemistryModel_H
#define LoadBalanced_CKJacChemistryModel_H

#include "clockTime.H"
#include "LoadBalancedChemistryModel.H"

// CKJac dependencies
extern "C" 
{
    #include "interface.H"
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class LoadBalanced_CKJacChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template <class ReactionThermo, class ThermoType>
class LoadBalanced_CKJacChemistryModel
    : public LoadBalancedChemistryModel<ReactionThermo, ThermoType> {
    template <class DeltaTType> scalar solve(const DeltaTType& deltaT);
		mutable scalar DTDt;
		mutable scalarList WDOT_;
		mutable scalarList CKJAC_;
public:
    //- Runtime type information
    TypeName("LB_CKJac");

    // Constructors

    //- Construct from thermo
    LoadBalanced_CKJacChemistryModel(const ReactionThermo& thermo);

    //- Destructor
    virtual ~LoadBalanced_CKJacChemistryModel();

    //- Disallow default bitwise copy construction
    LoadBalanced_CKJacChemistryModel(const LoadBalanced_CKJacChemistryModel&);

    inline virtual label nEqns() const;

    //- Override jacobian and derivatives function from StandardChemistryModel
    virtual void jacobian(const scalar        t,
                          const scalarField&  c,
                          const label         li,
                          scalarField&        dcdt,
                          scalarSquareMatrix& J) const override;

    virtual void derivatives(const scalar       t,
                             const scalarField& c,
                             const label        li,
                             scalarField&       dcdt) const override;

    //- Override Qdot function from StandardChemistryModel for using pyJac's enthalpy of formation
//    virtual tmp<volScalarField> Qdot() const override;

private:
    ///
    ///@brief List of enthalpy of formation from pyJac
    ///
    ///
//    scalarList sp_enth_form;

    //- Update the reaction rate of cell i
/*    void
    updateReactionRate(const ChemistrySolution& solution, const label& i);
*/    
    //- Return the main chemistry variable
    //  (concentration for standard, mass fraction for pyJac)
    scalarField getVariable(const scalarField& concentration, const scalarField& massFraction);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "LoadBalanced_CKJacChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "LoadBalanced_CKJacChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
