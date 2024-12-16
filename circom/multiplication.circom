template Multiplication() {
    signal input s1;
    signal input s2;
    signal output y;

    y <== s1 * s2;
}

component main = Multiplication();
