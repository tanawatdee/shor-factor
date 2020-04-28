using System;
using System.Threading.Tasks;

using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;

namespace Shor
{
    class Driver
    {
        static async Task Main(string[] args)
        {
            using (var qsim = new QuantumSimulator())
            {	// create each quantum circuit case
                int[,] arr = {
					{9, 4, 2}, //number to factor, no. of bits, random number
					{9, 4, 4},
					//...list of all possible cases
            	};
                int l = arr.GetLength(0);  
                // iterate through each cases
                for(int i = 0; i < l; i++){
                	ResourcesEstimator estimator = new ResourcesEstimator();
                	// run the machine class to estimate the resources used
	                await FindOrder.Run(estimator, arr[i,0], arr[i,1], arr[i,2]);
	                // query the data of depth, width, gates count
	                Console.WriteLine(
	                	"[" +
	                	arr[i,0] + "," +
	                	arr[i,2] + "," +
	                	estimator.Data.Rows.Find("CNOT")["Sum"] + "," +
	                	estimator.Data.Rows.Find("QubitClifford")["Sum"] + "," +
	                	estimator.Data.Rows.Find("R")["Sum"] + "," +
	                	estimator.Data.Rows.Find("Measure")["Sum"] + "," +
	                	estimator.Data.Rows.Find("T")["Sum"] + "," +
	                	estimator.Data.Rows.Find("Depth")["Sum"] + "," +
	                	estimator.Data.Rows.Find("Width")["Sum"] + "," +
	                	estimator.Data.Rows.Find("BorrowedWidth")["Sum"] + "],"
	                );
	            }
            }
        }
    }
}