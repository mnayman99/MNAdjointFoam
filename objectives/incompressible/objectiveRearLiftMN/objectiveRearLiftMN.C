/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019, 2022 PCOpt/NTUA
    Copyright (C) 2013-2019, 2022 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "objectiveRearLiftMN.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectiveRearLiftMN, 0);
addToRunTimeSelectionTable
(
    objectiveIncompressible,
    objectiveRearLiftMN,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectiveRearLiftMN::objectiveRearLiftMN
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveIncompressible(mesh, dict, adjointSolverName, primalSolverName),
    forcePatches_
    (
        mesh_.boundaryMesh().patchSet
        (
            dict.get<wordRes>("patches")
        ).sortedToc()
    ),
    forceDirection_(dict.get<vector>("liftDir")),
    momentDirection_(dict.get<vector>("axis")),
    rotationCentre_(dict.get<vector>("rotationCenter")),
    Aref_(dict.get<scalar>("Aref")),
    lRef_(dict.get<scalar>("lRef")),
    UInf_(dict.get<scalar>("UInf")),
    invDenom_(2./(UInf_*UInf_*Aref_)),
    invMomDenom_(2./(UInf_*UInf_*Aref_*lRef_)),
    devReff_(vars_.turbulence()->devReff()())
{
    // Sanity check and print info
    if (forcePatches_.empty())
    {
        FatalErrorInFunction
            << "No valid patch name on which to minimize " << type() << endl
            << exit(FatalError);
    }
    if (debug)
    {
        Info<< "Minimizing " << type() << " in patches:" << endl;
        for (const label patchI : forcePatches_)
        {
            Info<< "\t " << mesh_.boundary()[patchI].name() << endl;
        }
    }

    // Allocate boundary field pointers
    bdJdpPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdSdbMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdxdbMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdxdbDirectMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdJdnutPtr_.reset(createZeroBoundaryPtr<scalar>(mesh_));
    bdJdGradUPtr_.reset(createZeroBoundaryPtr<tensor>(mesh_));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar objectiveRearLiftMN::J()
{
    vector pressureForce(Zero);
    vector viscousForce(Zero);
    vector cumulativeForce(Zero);

    vector pressureMoment(Zero);
    vector viscousMoment(Zero);
    vector cumulativeMoment(Zero);


    const volScalarField& p = vars_.pInst();
    const autoPtr<incompressible::turbulenceModel>&
       turbulence = vars_.turbulence();

    volSymmTensorField devReff(turbulence->devReff());

    for (const label patchI : forcePatches_)
    {
        const vectorField& Sf = mesh_.Sf().boundaryField()[patchI];
        vectorField dx(mesh_.Cf().boundaryField()[patchI] - rotationCentre_);
        pressureForce += gSum(Sf*p.boundaryField()[patchI]);
        viscousForce += gSum(devReff.boundaryField()[patchI] & Sf);
        pressureMoment += gSum((dx ^ Sf)*p.boundaryField()[patchI]);
        viscousMoment += gSum((dx^(devReff_.boundaryField()[patchI] & Sf)));
    }

    cumulativeForce = pressureForce + viscousForce;
    cumulativeMoment = pressureMoment + viscousMoment;

    scalar force = cumulativeForce & forceDirection_;
    scalar moment = cumulativeMoment & momentDirection_;
    scalar rearLift = force - (moment/lRef_);

    // Intentionally not using denom - derived might implement virtual denom()
    // function differently
    scalar Cforce = force*invDenom_;
    scalar Cm = moment*invMomDenom_;
    scalar CLr = rearLift * invDenom_;

    DebugInfo
        << "Moment|Coeff " << moment << "|" << Cm << endl <<
        "Force|Coeff " << force << "|" << Cforce << endl <<
        "Rear Lift|Coeff " << rearLift << "|" << CLr << endl;

    J_ = CLr;

    return CLr;
}

void objectiveRearLiftMN::update_meanValues()
{
    if (computeMeanFields_)
    {
        const volVectorField& U = vars_.U();
        const autoPtr<incompressible::RASModelVariables>& turbVars =
            vars_.RASModelVariables();
        const singlePhaseTransportModel& lamTransp = vars_.laminarTransport();

        devReff_ = turbVars->devReff(lamTransp, U)();
    }
}

void objectiveRearLiftMN::update_boundarydJdp()
{
    for (const label patchI : forcePatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        tmp<vectorField> tdx = patch.Cf() - rotationCentre_;
        bdJdpPtr_()[patchI] = forceDirection_*invDenom_ - (momentDirection_ ^ tdx)*invMomDenom_;
    }
}


void objectiveRearLiftMN::update_dSdbMultiplier()
{
    // Compute contributions with mean fields, if present
    const volScalarField& p = vars_.p();
    
    for (const label patchI : forcePatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        tmp<vectorField> tdx = patch.Cf() - rotationCentre_;
        bdSdbMultPtr_()[patchI] =
        (
            (
                forceDirection_& devReff_.boundaryField()[patchI]
            )
          + (forceDirection_)*p.boundaryField()[patchI]
        )
       *invDenom_
       -
       (
            (
                (
                    (momentDirection_ ^ tdx()) &
                    (
                        devReff_.boundaryField()[patchI]
                    )
                )
            )
            + (momentDirection_ ^ tdx())*p.boundaryField()[patchI]
        )
        *invMomDenom_;
    }
}


void objectiveRearLiftMN::update_dxdbMultiplier()
{
    const volScalarField& p = vars_.p();
    const volVectorField& U = vars_.U();

    const autoPtr<incompressible::RASModelVariables>&
        turbVars = vars_.RASModelVariables();
    const singlePhaseTransportModel& lamTransp = vars_.laminarTransport();

    // We only need to modify the boundaryField of gradU locally.
    // If grad(U) is cached then
    // a. The .ref() call fails since the tmp is initialised from a
    //    const ref
    // b. we would be changing grad(U) for all other places in the code
    //    that need it
    // So, always allocate new memory and avoid registering the new field
    tmp<volTensorField> tgradU =
        volTensorField::New("gradULocal", fvc::grad(U));
    volTensorField& gradU = tgradU.ref();
    volTensorField::Boundary& gradUbf = gradU.boundaryFieldRef();

    // Explicitly correct the boundary gradient to get rid of
    // the tangential component
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        if (isA<wallFvPatch>(patch))
        {
            tmp<vectorField> tnf = patch.nf();
            gradUbf[patchI] = tnf*U.boundaryField()[patchI].snGrad();
        }
    }

    // Term coming from gradp
    tmp<volVectorField> tgradp(fvc::grad(p));
    const volVectorField& gradp = tgradp.cref();
    for (const label patchI : forcePatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        tmp<vectorField> tnf = patch.nf();
        tmp<vectorField> tdx = patch.Cf() - rotationCentre_;
        bdxdbMultPtr_()[patchI] =
            (forceDirection_ & mesh_.boundary()[patchI].nf())
           *gradp.boundaryField()[patchI]*invDenom_
           -
            (momentDirection_ & (tdx ^ tnf))*gradp.boundaryField()[patchI]
           *invMomDenom_;
    }
    tgradp.clear();

    // Term coming from stresses
    tmp<volScalarField> tnuEff = lamTransp.nu() + turbVars->nutRef();
    tmp<volSymmTensorField> tstress = tnuEff*twoSymm(tgradU);
    const volSymmTensorField& stress = tstress.cref();
    autoPtr<volVectorField> ptemp
        (Foam::createZeroFieldPtr<vector>( mesh_, "temp", sqr(dimVelocity)));
    volVectorField& temp = ptemp.ref();

    for (label idir = 0; idir < pTraits<vector>::nComponents; ++idir)
    {
        unzipRow(stress, idir, temp);
        volTensorField gradStressDir(fvc::grad(temp));
        for (const label patchI : forcePatches_)
        {
            const fvPatch& patch = mesh_.boundary()[patchI];
            tmp<vectorField> tnf = patch.nf();
            tmp<vectorField> tdx = patch.Cf() - rotationCentre_;
            tmp<scalarField> taux = (momentDirection_ ^ tdx)().component(idir);
            bdxdbMultPtr_()[patchI] -=
                forceDirection_.component(idir)
               *(gradStressDir.boundaryField()[patchI] & tnf)*invDenom_
               -
                taux*(gradStressDir.boundaryField()[patchI] & tnf)
               *invMomDenom_;
        }
    }
}


void objectiveRearLiftMN::update_dxdbDirectMultiplier()
{
    const volScalarField& p = vars_.p();

    for (const label patchI : forcePatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        tmp<vectorField> tnf = patch.nf();
        const vectorField& nf = tnf();
        const vectorField dx(patch.Cf() - rotationCentre_);
        const vectorField force
        (
            (
                ((p.boundaryField()[patchI]*nf)
              + (devReff_.boundaryField()[patchI] & nf))
            )
        );
        bdxdbDirectMultPtr_()[patchI] =
            -(force^momentDirection_)*invMomDenom_;
    }
}


void objectiveRearLiftMN::update_boundarydJdnut()
{
    const volVectorField& U = vars_.U();
    volSymmTensorField devGradU(dev(twoSymm(fvc::grad(U))));

    for (const label patchI : forcePatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        tmp<vectorField> tnf = patch.nf();
        tmp<vectorField> tdx = patch.Cf() - rotationCentre_;
        const fvPatchSymmTensorField& bdevGradU =
            devGradU.boundaryField()[patchI];
        bdJdnutPtr_()[patchI] =
          - ((bdevGradU & forceDirection_) & tnf)
           *invDenom_
           +
           ((tdx ^ (bdevGradU & tnf)) & momentDirection_)*invMomDenom_;
    }
}


void objectiveRearLiftMN::update_boundarydJdGradU()
{
    const autoPtr<incompressible::RASModelVariables>& turbVars =
        vars_.RASModelVariables();
    const singlePhaseTransportModel& lamTransp = vars_.laminarTransport();
    volScalarField nuEff(lamTransp.nu() + turbVars->nutRef());
    for (const label patchI : forcePatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        const vectorField& Sf = patch.Sf();
        bdJdGradUPtr_()[patchI] =
          - nuEff.boundaryField()[patchI]
           *dev(forceDirection_*Sf + Sf*forceDirection_);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //
