namespace Shor {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Bitwise;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Convert;

    operation FindOrder(N: Int, n: Int, Nr: Int) : Result[] {
        using (q = Qubit[4*n])  {
        	let q_sup = q[0*n..2*n-1];
		    let q_mod = q[2*n..3*n-1];
		    let q_acc = q[3*n..4*n-1];
		    //let q_anc = q[4*n..4*n-0];
            mutable Na = Nr;

		    ApplyToEach(H, q_sup); 	// uniform superposition
		    X(q_mod[0]);			// initial state |0..01>

		    for(i in 0..2*n-1){
                //out-of-place mod mul
                (Controlled MultiplyAndAddByModularInteger)([q_sup[i]], (Na, N, LittleEndian(q_mod), LittleEndian(q_acc)));
                //SWAP
                for(j in 0..n-1){
                    SWAP(q_mod[j], q_acc[j]);
                }
                //mod inverse and reverse
                (Adjoint Controlled MultiplyAndAddByModularInteger)([q_sup[i]], (InverseModI(Na, N), N, LittleEndian(q_mod), LittleEndian(q_acc)));
                //Message(IntAsStringWithFormat(Na, "{0}"));
                set Na = ModulusI(PowI(Na, 2), N);
            }

            QFTLE(LittleEndian(q_sup));

            let r = MultiM(q_sup);
            ResetAll(q);
            return r;
        }
    }
}