pub trait GPUScalarMul {

}

pub trait GPUScalarMulSlice {

}

impl<G: AffineCurve> GPUScalarMulSlice for [G] {
    
}

pub trait GPUScalarMulParameters {

}
