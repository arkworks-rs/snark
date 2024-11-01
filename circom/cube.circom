template Cube() {
    
    signal input x;
    signal x2;
    
    x2 <== x * x;

    x * x2 === 27;
}

component main = Cube();
