    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhiH, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
//        sources(rho, U)
        fvOptions(rho, U)
    );

    UEqn.relax();
    fvOptions.constrain(UEqn);

//    sources.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
//			mixture.surfaceTensionForce()
                    sigma*fvc::snGrad(psi)*fvc::interpolate(C)*fvc::interpolate(delta)
		  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
	 fvOptions.correct(U);
    }
