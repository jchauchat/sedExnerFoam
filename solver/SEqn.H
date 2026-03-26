{
    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + turbulence->nut()/Prt);

    fvScalarMatrix SEqn
    (
        fvm::ddt(S)
      + fvm::div(phi, S)
      - fvm::laplacian(alphaEff, S)
     ==
        fvOptions(S)
    );

    SEqn.relax();

    fvOptions.constrain(SEqn);

    SEqn.solve();

    fvOptions.correct(S);

    // Add Boussinesq term
    rhok = 1.0 + Cs * (rhoS - rhoF) / rhoF + betaS * S;
}
