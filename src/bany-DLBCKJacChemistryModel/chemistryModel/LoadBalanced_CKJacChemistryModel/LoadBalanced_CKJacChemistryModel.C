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

\*---------------------------------------------------------------------------*/

#include "LoadBalanced_CKJacChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam {

template <class ReactionThermo, class ThermoType>
LoadBalanced_CKJacChemistryModel<ReactionThermo, ThermoType>::LoadBalanced_CKJacChemistryModel(
    const ReactionThermo& thermo)
    : LoadBalancedChemistryModel<ReactionThermo, ThermoType>(thermo), 
//	sp_enth_form(this->nSpecie_),
	DTDt(0),
	WDOT_(this->nSpecie_),
	CKJAC_((this->nSpecie_+1)*(this->nSpecie_+1))
	{
    
    Info << "Overriding StandardChemistryModel by LoadBalanced_CKJacChemistryModel:" << endl; 

	Info << "CKJac mechanism information:" << 
                "\n\tNumber of species: " << this->nSpecie_ <<
                "\n\tNumber of forward reactions: " << this->nReaction_ << "\n" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class ReactionThermo, class ThermoType>
LoadBalanced_CKJacChemistryModel<ReactionThermo, ThermoType>::~LoadBalanced_CKJacChemistryModel() {}

template<class ReactionThermo, class ThermoType>
void LoadBalanced_CKJacChemistryModel<ReactionThermo, ThermoType>::jacobian
(
	const scalar t, 
	const scalarField& c, 
	const label li, 
	scalarField& dcdt, 
	scalarSquareMatrix& J
) const 
{
//this->numbofJac++;
    const scalar T = c[this->nSpecie_];
    const scalar p = c[this->nSpecie_ + 1];

    forAll(this->c_, i)
    {
        this->c_[i] = max(c[i], 0);
    }
	
    J = Zero;
    dcdt = Zero;

    forAll(this->c_, i)
    {
        this->c_[i] = 1E-3*max(c[i], 0); // kmol/m3 -> mol/cm3
//		Info << "this->c_["<<i<<"] = " << this->c_[i] << endl;
    }

	CALL_JAC_OF_C(&this->nSpecie_, &p, &T, &this->c_[0], &CKJAC_[0]);

    forAll(this->c_, i)
    {
        this->c_[i] *= 1E3; // mol/cm3 -> kmol/m3
    }

// J[i][j]        = PARTIAL f(i)/PARTIAL C(j)
// J[i][nSpecie_] = PARTIAL f(i)/PARTIAL temp
// J[nSpecie_][j] = PARTIAL f(temp)/PARTIAL C(j)
// J[nSpecie_][nSpecie_] = PARTIAL f(temp)/PARTIAL temp
// J[nSpecie_+1][j] = PARTIAL f(pressure)/PARTIAL C(j) = 0
// J[nSpecie_+2][j] = PARTIAL f(discharge)/PARTIAL C(j)
// J[nSpecie_+2][nSpecie_] = PARTIAL f(discharge)/PARTIAL temp
	
// PARTIAL f(Ci)/PARTIAL f(Cj)
    for (label i = 0; i < this->nSpecie_; i++)
    {
		for (label j = 0; j < this->nSpecie_; j++)
		{
			J[i][j] = CKJAC_[(this->nSpecie_+1)*(i+1)+(j+1)];  // 1/s
		}
    }

// PARTIAL f(Cj)/PARTIAL T
	for (label j = 0; j < this->nSpecie_; j++)
	{
		J[j][this->nSpecie_] = CKJAC_[(this->nSpecie_+1)*(j+1)] * 1E3; // kmol/(m3*s*K)
	}

// PARTIAL f(temp)/PARTIAL C(j)
	for (label j = 0; j < this->nSpecie_; j++)
	{
		J[this->nSpecie_][j] = CKJAC_[j+1] * 1E-3; // K*m3/(kmol*s)
	}

// PARTIAL f(temp)/PARTIAL T 
	J[this->nSpecie_][this->nSpecie_] = CKJAC_[0]; // 1/s
	
/*	
	label numbofzero = 0;
	for (label i = 0; i < this->nSpecie_+2; i++)
    {
		for (label j = 0; j < this->nSpecie_+2; j++)
		{
			if (J(i,j)==0)
			{
				numbofzero++;
			}
		}
    } 
	this->sparse = numbofzero*0.1/((this->nSpecie_+2)*(this->nSpecie_+2)*0.1);
*/	
/*	Info <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	for (label i = 0; i < this->nSpecie_+2; i++)
    {
		Info << i;
		for (label j = 0; j < this->nSpecie_+2; j++)
		{
			Info <<" "<< J(i, j)<<"  ";
		}
		Info <<endl;
    } 
	Info <<"---------------------------------------------------------------"<<endl;
*/

}

template<class ReactionThermo, class ThermoType>
void LoadBalanced_CKJacChemistryModel<ReactionThermo, ThermoType>::derivatives
(
	const scalar t, 
	const scalarField& c, 
	const label li, 
	scalarField& dcdt
) const 
{

//	Info << setprecision(15);
//	Info <<  "in CKChemistryModel"  << endl;
    const scalar T = c[this->nSpecie_];
    const scalar p = c[this->nSpecie_ + 1];
//	Info << "T = " << T << endl;
//	Info << "p = " << p << endl;

//Info << "c = " << c << endl;

    forAll(this->c_, i)
    {
        this->c_[i] = 1E-3*max(c[i], 0); // kmol/m3 -> mol/cm3
    }
	
	CALL_OMEGA_OF_C(&T, &this->c_[0], &p, &this->nSpecie_, &WDOT_[0], &DTDt); // WDOT_ mol/(cm3*s)
	
	for (label i = 0; i < this->nSpecie_; i++)
	{
		dcdt[i] = WDOT_[i]*1E3; // mol/(cm3*s) -> kmol/(m3*s) // The specie index in OF and fortran must be consistent !!!
		this->c_[i] *= 1E3; // mol/cm3 -> kmol/m3
//		Info << "dcdt["<<i<<"] = " << dcdt[i] << endl;
	}

	// dT/dt
	dcdt[this->nSpecie_] = DTDt;
//Info << "DTDt = " << DTDt << endl;
    // dp/dt = ...
    dcdt[this->nSpecie_ + 1] = 0;
	
}

/*
// TODO: Make this work for pyjac
template <class ReactionThermo, class ThermoType>
void Foam::LoadBalanced_CKJacChemistryModel<ReactionThermo, ThermoType>::updateReactionRate
(
    const ChemistrySolution& solution, const label& i
)
{
    for(label j = 0; j < this->nSpecie_; j++)
    {
        this->RR_[j][i] = solution.c_increment[j] * solution.rhoi;
    }
    this->deltaTChem_[i] = min(solution.deltaTChem, this->deltaTChemMax_);
}
*/
template <class ReactionThermo, class ThermoType>
Foam::scalarField Foam::LoadBalanced_CKJacChemistryModel<ReactionThermo, ThermoType>::getVariable
(
    const scalarField& concentration, const scalarField& massFraction
)
{
    return concentration;
}


} // namespace Foam
