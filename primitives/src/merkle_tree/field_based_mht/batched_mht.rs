/*
impl<F: PrimeField + MulShort, P: PoseidonParameters<Fr = F>> BatchFieldBasedHash for PoseidonBatchHash<F, P> {
    type Data = F;
    type Parameters = P;

    fn batch_evaluate(input_array: &mut[F], output_array: &mut[F]) {

        // Input:
        // This function calculates the hashes of inputs by groups of the rate P::R.
        // The inputs are arranged in an array and arranged as consecutive groups
        // Example:
        // For P::R = 2,
        // (d_00, d01, d_10, d_11, d_20, d_21, ...
        // Output:
        // The output will be placed in the output_array taking input length / P::R

        // Checks that size of input/output vector
        let array_length = input_array.len() / P::R;
        assert_eq!(output_array.len() >= array_length, true, "Not enough space for output vector compared to the input vector.");

        // Assign pre-computed values of the state vector equivalent to a permutation with zero element state vector
        //let state_z = vec![P::AFTER_ZERO_PERM[0], P::AFTER_ZERO_PERM[1], P::AFTER_ZERO_PERM[2]];
        let mut state_z = Vec::new();
        for i in 0..P::T {
            state_z.push(P::AFTER_ZERO_PERM[i]);
        }

        // Copy the result of the permutation to a vector of state vectors of the length equal to the length of the input
        // state is a vector of 3-element state vector.
        let mut state = Vec::new();
        for _i in 0..array_length {
            state.push(state_z.clone());
        }

        // input_idx is to scan the input_array
        let mut input_idx = 0;

        for k in 0..array_length {
            for j in 0..P::R {
                state[k][j] += &input_array[input_idx];
                input_idx += 1;
            }
            // constant for m-ary Merkle tree
            state[k][P::R] += &P::C2;
        }

        // apply permutation after adding the input vector
        Self::poseidon_perm_gen(&mut state);

        // write the result of the hash extracted from the state vector to the output vector
        for k in 0..array_length {
            output_array[k] = state[k][0];
        }
    }

    fn merkle_tree(input_vec: &mut[Self::Data], output_vec: &mut[Self::Data], input_size: usize){
        // Supporting function that processes the inputs and outputs in chunks

        let num_cores = 16;

        assert_eq!(input_vec.len() % 2, 0, "The length of the input to the hash is not even.");
        assert_eq!(output_vec.len() >= input_vec.len() / 2, true,  "The length of the output is not greater or equal to half of the input length.");

        if input_size < 2 * num_cores {
            input_vec.par_chunks_mut(2).zip(output_vec.par_chunks_mut(1)).for_each( |(p1,p2)| {
                Self::batch_evaluate(p1, p2);
            });
            return;
        }

        use rayon::prelude::*;

        input_vec.par_chunks_mut(input_size/num_cores).zip(output_vec.par_chunks_mut(input_size/num_cores/2)).for_each( |(p1, p2)| {
            Self::batch_evaluate(p1, p2);
        });

    }

    fn merkle_tree_2_1(array: &mut [Self::Data], input_size: usize) {
        // Main function for the Merkle-tree computation

        let mut copy_vec = &mut array[..];
        let mut size_input = input_size;

        while size_input > 1{
            let (input_vec, output_vec) = copy_vec.split_at_mut(size_input);
            Self::merkle_tree(input_vec, output_vec, size_input);
            copy_vec = output_vec;
            size_input = size_input / 2;
        }
    }

}*/
