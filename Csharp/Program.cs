 using System;
 using System.IO;

/* Program By Amin Khiabani
   aminkhiabani@outlook.com
*/

namespace LegendreFunction
{
    class Program
    {
        static void Main(string[] args)
        {
            double normal_g = 0.0; 
            double adopted_g = 0.0;
            double calc_g0 = 0.0;
            double calc_g1 = 0.0;

            Console.WriteLine();
            Console.WriteLine("           normal        Adopted         Calculated        Calculated ");
            Console.WriteLine("phi        gravity       Gravity          gravity0          gravity1 ");
            for (double ph=0; ph<=90; ph = ph + 5.0)
            {
                
                normal_g =  Legendre.normalGravity(ph)*1e5;
                adopted_g = Legendre.adaptedGravity(ph, 0.0)*1e5;
                calc_g0 = Legendre.Calcgravity(ph, 0)*1e5;
                calc_g1 = Legendre.Calcgravity(ph, 1)*1e5;
                Console.WriteLine("{0:0.0}    {1:0.000000}    {2:0.000000}    {3:0.000000}    {4:0.000000}",
                    ph,normal_g, adopted_g,calc_g0,calc_g1);                
            }
            
         /*   
            string fileName = @"file directory\dV_ELL_Earth2014.gfc";
                if (File.Exists(fileName))
                {
                    double Phi = 45.0;
                    double[] Lambda = { 46.0, 47.0, 48.0 };
                    for(int i=0; i<=2; i++)
                    {
                        Console.WriteLine();
                        double del_g = Legendre.gravitation(Phi, Lambda[i], 0, fileName);                    
                        double calc_g = Legendre.Calcgravity(Phi, 1);
                        Console.WriteLine(" Delta_Gravity at Latitude {0} , Longitude {1} is {2}",
                            Phi, Lambda[i], del_g);
                        double g = (calc_g + del_g) * 1E5; 
                        Console.WriteLine(" Calculated Gravity at Latitude {0} , Longitude {1} is {2}",
                            Phi, Lambda[i], g);
                    }
               
                }
                else
                {
                    Console.WriteLine("File Does Not Exists");
                }
          */  
        }                                
    }
}
