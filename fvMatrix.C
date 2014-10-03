template<class Type>
void Foam::fvMatrix<Type>::relax2(const scalar alpha)
{
    if (alpha <= 0)
    {
        return;
    }

    if (debug)
    {
        InfoIn("fvMatrix<Type>::relax(const scalar alpha)")
            << "Relaxing " << psi_.name() << " by " << alpha
            << endl;
    }

    Field<Type>& S = source();
    scalarField& D = diag();

    // Store the current unrelaxed diagonal for use in updating the source
    scalarField D0(D);

    // Calculate the sum-mag off-diagonal from the interior faces
    scalarField sumOff(D.size(), 0.0);
    sumMagOffDiag(sumOff);

    // Handle the boundary contributions to the diagonal
    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.size())
        {
            const labelUList& pa = lduAddr().patchAddr(patchI);
            Field<Type>& iCoeffs = internalCoeffs_[patchI];

            if (ptf.coupled())
            {
                const Field<Type>& pCoeffs = boundaryCoeffs_[patchI];

                // For coupled boundaries add the diagonal and
                // off-diagonal contributions
                forAll(pa, face)
                {
                    D[pa[face]] += component(iCoeffs[face], 0);
                    sumOff[pa[face]] += mag(component(pCoeffs[face], 0));
                }
            }
            else
            {
                // For non-coupled boundaries add the maximum magnitude diagonal
                // contribution to ensure stability
                forAll(pa, face)
                {
                    D[pa[face]] += cmptMax(cmptMag(iCoeffs[face]));
                }
            }
        }
    }


    if (debug)
    {
        // Calculate amount of non-dominance.
        label nNon = 0;
        scalar maxNon = 0.0;
        scalar sumNon = 0.0;
        forAll(D, celli)
        {
            scalar d = (sumOff[celli] - D[celli])/mag(D[celli]);

            if (d > 0)
            {
                nNon++;
                maxNon = max(maxNon, d);
                sumNon += d;
            }
        }

        reduce(nNon, sumOp<label>());
        reduce(maxNon, maxOp<scalar>());
        reduce(sumNon, sumOp<scalar>());
        sumNon /= returnReduce(D.size(), sumOp<label>());

        InfoIn("fvMatrix<Type>::relax(const scalar alpha)")
            << "Matrix dominance test for " << psi_.name() << nl
            << "    number of non-dominant cells   : " << nNon << nl
            << "    maximum relative non-dominance : " << maxNon << nl
            << "    average relative non-dominance : " << sumNon << nl
            << endl;
    }


    // Ensure the matrix is diagonally dominant...
    // Assumes that the central coefficient is positive and ensures it is
    forAll(D, celli)
    {
        D[celli] = max(mag(D[celli]), sumOff[celli]);
    }

    // ... then relax
    D /= alpha;

    // Now remove the diagonal contribution from coupled boundaries
    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.size())
        {
            const labelUList& pa = lduAddr().patchAddr(patchI);
            Field<Type>& iCoeffs = internalCoeffs_[patchI];

            if (ptf.coupled())
            {
                forAll(pa, face)
                {
                    D[pa[face]] -= component(iCoeffs[face], 0);
                }
            }
            else
            {
                forAll(pa, face)
                {
                    D[pa[face]] -= cmptMin(iCoeffs[face]);
                }
            }
        }
    }

    // Finally add the relaxation contribution to the source.
    S += (D - D0)*psi_.internalField();
}
