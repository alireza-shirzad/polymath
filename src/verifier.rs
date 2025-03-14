use ark_ec::{pairing::Pairing, AffineRepr, VariableBaseMSM};
use ark_ff::PrimeField;
use ark_std::One;

use crate::{
    common::{B_POLYMATH, MINUS_ALPHA, MINUS_GAMMA},
    Polymath, PolymathError, Transcript, VerifyingKey,
};

use super::Proof;

impl<F: PrimeField, E, T> Polymath<E, T>
where
    E: Pairing<ScalarField = F>,
    T: Transcript<Challenge = F>,
{
    /// Verify a Polymath proof `proof` against the verification key `vk`,
    /// with respect to the instance `public_inputs`.
    pub(crate) fn verify_proof(
        vk: &VerifyingKey<E>,
        proof: &Proof<E>,
        public_inputs: &[F],
    ) -> Result<bool, PolymathError> {
        let mut t = T::new(B_POLYMATH);

        let public_inputs = &[&[F::one()], public_inputs].concat();

        // compute challenge x1
        let x1: F = Self::compute_x1(&mut t, public_inputs, &[proof.a_g1, proof.c_g1])?;

        // compute y1=x1^sigma
        let y1: F = Self::compute_y1(x1, vk.sigma);

        let y_inverse = y1.inverse().unwrap();
        let y1_gamma = y_inverse.pow([MINUS_GAMMA]);
        let y1_alpha = y_inverse.pow([MINUS_ALPHA]);
        let pi_at_x1 = Self::compute_pi_at_x1(vk, public_inputs, x1, y1_gamma);

        // compute c_at_x1
        let c_at_x1 = Self::compute_c_at_x1(y1_gamma, y1_alpha, proof.a_at_x1, pi_at_x1);

        let x2 = Self::compute_x2(&mut t, &x1, &[proof.a_at_x1, c_at_x1])?;

        let commitments_minus_evals_in_g1 = proof.a_g1 - proof.c_g1 * x2 - vk.e.one_g1 * (proof.a_at_x1 + x2 * c_at_x1);
        
        let x_minus_x1_in_g2 = vk.e.x_g2 - vk.e.one_g2 * x1;
        
        let pairing_output = E::multi_pairing(
            [
                commitments_minus_evals_in_g1.into(),
                (-(proof.d_g1.into_group())).into(),
            ],
            [
                vk.e.z_g2.clone(),
                x_minus_x1_in_g2.into().into(),
            ],
        );

        Ok(pairing_output.0.is_one())
    }
}
