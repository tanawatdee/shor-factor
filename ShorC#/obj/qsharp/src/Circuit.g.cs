//------------------------------------------------------------------------------
// <auto-generated>                                                             
//     This code was generated by a tool.                                       
//     Changes to this file may cause incorrect behavior and will be lost if    
//     the code is regenerated.                                                 
// </auto-generated>                                                            
//------------------------------------------------------------------------------
#pragma warning disable 162
#pragma warning disable 1591
using System;
using Microsoft.Quantum.Core;
using Microsoft.Quantum.Intrinsic;
using Microsoft.Quantum.Simulation.Core;

[assembly: CallableDeclaration("{\"Kind\":{\"Case\":\"Operation\"},\"QualifiedName\":{\"Namespace\":\"Shor\",\"Name\":\"FindOrder\"},\"Attributes\":[],\"SourceFile\":\"F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs\",\"Position\":{\"Item1\":9,\"Item2\":4},\"SymbolRange\":{\"Item1\":{\"Line\":1,\"Column\":11},\"Item2\":{\"Line\":1,\"Column\":20}},\"ArgumentTuple\":{\"Case\":\"QsTuple\",\"Fields\":[[{\"Case\":\"QsTupleItem\",\"Fields\":[{\"VariableName\":{\"Case\":\"ValidName\",\"Fields\":[\"N\"]},\"Type\":{\"Case\":\"Int\"},\"InferredInformation\":{\"IsMutable\":false,\"HasLocalQuantumDependency\":false},\"Position\":{\"Case\":\"Null\"},\"Range\":{\"Item1\":{\"Line\":1,\"Column\":21},\"Item2\":{\"Line\":1,\"Column\":22}}}]},{\"Case\":\"QsTupleItem\",\"Fields\":[{\"VariableName\":{\"Case\":\"ValidName\",\"Fields\":[\"n\"]},\"Type\":{\"Case\":\"Int\"},\"InferredInformation\":{\"IsMutable\":false,\"HasLocalQuantumDependency\":false},\"Position\":{\"Case\":\"Null\"},\"Range\":{\"Item1\":{\"Line\":1,\"Column\":29},\"Item2\":{\"Line\":1,\"Column\":30}}}]},{\"Case\":\"QsTupleItem\",\"Fields\":[{\"VariableName\":{\"Case\":\"ValidName\",\"Fields\":[\"Nr\"]},\"Type\":{\"Case\":\"Int\"},\"InferredInformation\":{\"IsMutable\":false,\"HasLocalQuantumDependency\":false},\"Position\":{\"Case\":\"Null\"},\"Range\":{\"Item1\":{\"Line\":1,\"Column\":37},\"Item2\":{\"Line\":1,\"Column\":39}}}]}]]},\"Signature\":{\"TypeParameters\":[],\"ArgumentType\":{\"Case\":\"TupleType\",\"Fields\":[[{\"Case\":\"Int\"},{\"Case\":\"Int\"},{\"Case\":\"Int\"}]]},\"ReturnType\":{\"Case\":\"ArrayType\",\"Fields\":[{\"Case\":\"Result\"}]},\"Information\":{\"Characteristics\":{\"Case\":\"EmptySet\"},\"InferredInformation\":{\"IsSelfAdjoint\":false,\"IsIntrinsic\":false}}},\"Documentation\":[]}")]
[assembly: SpecializationDeclaration("{\"Kind\":{\"Case\":\"QsBody\"},\"TypeArguments\":{\"Case\":\"Null\"},\"Information\":{\"Characteristics\":{\"Case\":\"EmptySet\"},\"InferredInformation\":{\"IsSelfAdjoint\":false,\"IsIntrinsic\":false}},\"Parent\":{\"Namespace\":\"Shor\",\"Name\":\"FindOrder\"},\"Attributes\":[],\"SourceFile\":\"F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs\",\"Position\":{\"Item1\":9,\"Item2\":4},\"HeaderRange\":{\"Item1\":{\"Line\":1,\"Column\":11},\"Item2\":{\"Line\":1,\"Column\":20}},\"Documentation\":[]}")]
#line hidden
namespace Shor
{
    [SourceLocation("F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs", OperationFunctor.Body, 10, -1)]
    public partial class FindOrder : Operation<(Int64,Int64,Int64), IQArray<Result>>, ICallable
    {
        public FindOrder(IOperationFactory m) : base(m)
        {
        }

        public class In : QTuple<(Int64,Int64,Int64)>, IApplyData
        {
            public In((Int64,Int64,Int64) data) : base(data)
            {
            }

            System.Collections.Generic.IEnumerable<Qubit> IApplyData.Qubits => null;
        }

        String ICallable.Name => "FindOrder";
        String ICallable.FullName => "Shor.FindOrder";
        public static OperationInfo<(Int64,Int64,Int64), IQArray<Result>> Info => new OperationInfo<(Int64,Int64,Int64), IQArray<Result>>(typeof(FindOrder));
        protected ICallable<IQArray<Qubit>, Microsoft.Quantum.Arithmetic.LittleEndian> MicrosoftQuantumArithmeticLittleEndian
        {
            get;
            set;
        }

        protected IUnitary<(Int64,Int64,Microsoft.Quantum.Arithmetic.LittleEndian,Microsoft.Quantum.Arithmetic.LittleEndian)> MicrosoftQuantumArithmeticMultiplyAndAddByModularInteger
        {
            get;
            set;
        }

        protected ICallable MicrosoftQuantumCanonApplyToEach
        {
            get;
            set;
        }

        protected IUnitary<Microsoft.Quantum.Arithmetic.LittleEndian> MicrosoftQuantumCanonQFTLE
        {
            get;
            set;
        }

        protected Allocate Allocate
        {
            get;
            set;
        }

        protected IUnitary<Qubit> MicrosoftQuantumIntrinsicH
        {
            get;
            set;
        }

        protected Release Release
        {
            get;
            set;
        }

        protected ICallable<IQArray<Qubit>, QVoid> MicrosoftQuantumIntrinsicResetAll
        {
            get;
            set;
        }

        protected IUnitary<(Qubit,Qubit)> MicrosoftQuantumIntrinsicSWAP
        {
            get;
            set;
        }

        protected IUnitary<Qubit> MicrosoftQuantumIntrinsicX
        {
            get;
            set;
        }

        protected ICallable<(Int64,Int64), Int64> MicrosoftQuantumMathInverseModI
        {
            get;
            set;
        }

        protected ICallable<(Int64,Int64), Int64> MicrosoftQuantumMathModulusI
        {
            get;
            set;
        }

        protected ICallable<(Int64,Int64), Int64> MicrosoftQuantumMathPowI
        {
            get;
            set;
        }

        protected ICallable<IQArray<Qubit>, IQArray<Result>> MicrosoftQuantumMeasurementMultiM
        {
            get;
            set;
        }

        public override Func<(Int64,Int64,Int64), IQArray<Result>> Body => (__in__) =>
        {
            var (N,n,Nr) = __in__;
#line hidden
            {
#line 11 "F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs"
                var q = Allocate.Apply((4L * n));
#line hidden
                bool __arg1__ = true;
                try
                {
#line 12 "F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs"
                    var q_sup = (IQArray<Qubit>)q?.Slice(new QRange((0L * n), ((2L * n) - 1L)));
#line 13 "F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs"
                    var q_mod = (IQArray<Qubit>)q?.Slice(new QRange((2L * n), ((3L * n) - 1L)));
#line 14 "F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs"
                    var q_acc = (IQArray<Qubit>)q?.Slice(new QRange((3L * n), ((4L * n) - 1L)));
#line 16 "F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs"
                    var Na = Nr;
#line 18 "F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs"
                    MicrosoftQuantumCanonApplyToEach.Apply((MicrosoftQuantumIntrinsicH, q_sup));
#line 19 "F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs"
                    MicrosoftQuantumIntrinsicX.Apply(q_mod[0L]);
#line 21 "F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs"
                    foreach (var i in new QRange(0L, ((2L * n) - 1L)))
#line hidden
                    {
#line 23 "F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs"
                        MicrosoftQuantumArithmeticMultiplyAndAddByModularInteger.Controlled.Apply((new QArray<Qubit>(q_sup[i]), (Na, N, new Microsoft.Quantum.Arithmetic.LittleEndian(q_mod), new Microsoft.Quantum.Arithmetic.LittleEndian(q_acc))));
#line 25 "F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs"
                        foreach (var j in new QRange(0L, (n - 1L)))
#line hidden
                        {
#line 26 "F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs"
                            MicrosoftQuantumIntrinsicSWAP.Apply((q_mod[j], q_acc[j]));
                        }

#line 29 "F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs"
                        MicrosoftQuantumArithmeticMultiplyAndAddByModularInteger.Controlled.Adjoint.Apply((new QArray<Qubit>(q_sup[i]), (MicrosoftQuantumMathInverseModI.Apply((Na, N)), N, new Microsoft.Quantum.Arithmetic.LittleEndian(q_mod), new Microsoft.Quantum.Arithmetic.LittleEndian(q_acc))));
#line 31 "F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs"
                        Na = MicrosoftQuantumMathModulusI.Apply((MicrosoftQuantumMathPowI.Apply((Na, 2L)), N));
                    }

#line 34 "F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs"
                    MicrosoftQuantumCanonQFTLE.Apply(new Microsoft.Quantum.Arithmetic.LittleEndian(q_sup));
#line 36 "F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs"
                    var r = (IQArray<Result>)MicrosoftQuantumMeasurementMultiM.Apply(q_sup);
#line 37 "F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs"
                    MicrosoftQuantumIntrinsicResetAll.Apply(q);
#line 38 "F:/OneDrive/Y4/Q%23/ShorDepth/Circuit.qs"
                    return r;
                }
#line hidden
                catch
                {
                    __arg1__ = false;
                    throw;
                }
#line hidden
                finally
                {
                    if (__arg1__)
                    {
#line hidden
                        Release.Apply(q);
                    }
                }
            }
        }

        ;
        public override void Init()
        {
            this.MicrosoftQuantumArithmeticLittleEndian = this.Factory.Get<ICallable<IQArray<Qubit>, Microsoft.Quantum.Arithmetic.LittleEndian>>(typeof(Microsoft.Quantum.Arithmetic.LittleEndian));
            this.MicrosoftQuantumArithmeticMultiplyAndAddByModularInteger = this.Factory.Get<IUnitary<(Int64,Int64,Microsoft.Quantum.Arithmetic.LittleEndian,Microsoft.Quantum.Arithmetic.LittleEndian)>>(typeof(Microsoft.Quantum.Arithmetic.MultiplyAndAddByModularInteger));
            this.MicrosoftQuantumCanonApplyToEach = this.Factory.Get<ICallable>(typeof(Microsoft.Quantum.Canon.ApplyToEach<>));
            this.MicrosoftQuantumCanonQFTLE = this.Factory.Get<IUnitary<Microsoft.Quantum.Arithmetic.LittleEndian>>(typeof(Microsoft.Quantum.Canon.QFTLE));
            this.Allocate = this.Factory.Get<Allocate>(typeof(Microsoft.Quantum.Intrinsic.Allocate));
            this.MicrosoftQuantumIntrinsicH = this.Factory.Get<IUnitary<Qubit>>(typeof(Microsoft.Quantum.Intrinsic.H));
            this.Release = this.Factory.Get<Release>(typeof(Microsoft.Quantum.Intrinsic.Release));
            this.MicrosoftQuantumIntrinsicResetAll = this.Factory.Get<ICallable<IQArray<Qubit>, QVoid>>(typeof(Microsoft.Quantum.Intrinsic.ResetAll));
            this.MicrosoftQuantumIntrinsicSWAP = this.Factory.Get<IUnitary<(Qubit,Qubit)>>(typeof(Microsoft.Quantum.Intrinsic.SWAP));
            this.MicrosoftQuantumIntrinsicX = this.Factory.Get<IUnitary<Qubit>>(typeof(Microsoft.Quantum.Intrinsic.X));
            this.MicrosoftQuantumMathInverseModI = this.Factory.Get<ICallable<(Int64,Int64), Int64>>(typeof(Microsoft.Quantum.Math.InverseModI));
            this.MicrosoftQuantumMathModulusI = this.Factory.Get<ICallable<(Int64,Int64), Int64>>(typeof(Microsoft.Quantum.Math.ModulusI));
            this.MicrosoftQuantumMathPowI = this.Factory.Get<ICallable<(Int64,Int64), Int64>>(typeof(Microsoft.Quantum.Math.PowI));
            this.MicrosoftQuantumMeasurementMultiM = this.Factory.Get<ICallable<IQArray<Qubit>, IQArray<Result>>>(typeof(Microsoft.Quantum.Measurement.MultiM));
        }

        public override IApplyData __dataIn((Int64,Int64,Int64) data) => new In(data);
        public override IApplyData __dataOut(IQArray<Result> data) => data;
        public static System.Threading.Tasks.Task<IQArray<Result>> Run(IOperationFactory __m__, Int64 N, Int64 n, Int64 Nr)
        {
            return __m__.Run<FindOrder, (Int64,Int64,Int64), IQArray<Result>>((N, n, Nr));
        }
    }
}