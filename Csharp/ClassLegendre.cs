using System;
using System.Collections.Generic;
using System.Text;

/* Program By Amin Khiabani
   aminkhiabani@outlook.com
*/

namespace LegendreFunction
{
    static class Legendre
    {
        public static double PI = 3.1415926535897932384626430;
        public static double DEGRAD = PI / 180.0;
        public static double a = 6378137.0;
        public static double ep2 = 0.00673949677540;
        public static double Omega = 0.7292115E-4;
        public static double GM = 0.39860050000E+15;
        public static double e2 = 0.006694380022900;
        public static double k = 0.0019318513530;

        public static int factorial(int n)
        {
            int f = 1;
            if (n == 0)
            {
                return f;
            }
            else
            {
                f = n * factorial(n - 1);
            }
            return f;
        }

        public static double LegendreF(int n, double x)
        {
            double LFS0 = 1.0;
            double LFS1 = x;
            double Lf = 0.0;
            if(n == 0)
            {
                Lf = LFS0;
            }
            else if(n == 1)
            {
               Lf = LFS1;
            }
            else if (n >= 1)
            {
                for (int i = 2; i <=n; i++)
                {
                    Lf = (-(i - 1) * LFS0 + (2 * i - 1) * x * LFS1) / i;
                    LFS0 = LFS1;
                    LFS1 = Lf;
                }                   
            }
            return Lf;
        }

        public static double LegendreR(int n, double x)
        {
            double Pn = 0.0;
        
         if (n == 0)
            {
                return 1.0;
            }           
         else if (n == 1)
            {
                return x;
            }           
         else
            {
                Pn = (-(n-1)*LegendreR(n-2,x) + (2*n-1)*x*LegendreR(n-1,x))/n;
            }

            return Pn;                 
        }

        public static double[] LegendreSeries(int n , double x)
        {
            double[] LFS = new double[n+1];

            LFS[0] = 1.0;
            LFS[1] = x;
            for (int I = 2; I <= n; I++)
                {
                    LFS[I] = (-(I - 1) * LFS[I - 2] + (2 * I - 1) * x * LFS[I - 1]) / I;
                }                        
            return LFS;
        }

        public static double[] ALegendFnmSeri(int n, double x)
        {
            double[] Pnm = new double[n+1];
            double[] ALFS = new double[n+1];
            double[] BLFS = new double[n+1];            

            double s = Math.Sqrt (1.0 - x * x);
           
            switch (n)
            {
                case 0 :
                    Pnm[0] = 1.0;                    
                    break;
                case 1 :
                    Pnm[0] = x;
                    Pnm[1] = s;
                    break;
                case 2 :
                    Pnm[0] = (3.0*x*x-1.0)/2.0;
                    Pnm[1] = 3.0 * s * x;
                    Pnm[2] = 3.0 * (1.0 - x * x);
                    break;
                case 3 :
                    Pnm[0] = x * (5.0 * x * x - 3.0) / 2.0;
                    Pnm[1] = 3.0 * (5.0 * x * x - 1.0) * s / 2.0;
                    Pnm[2] = 15.0 * x * (1.0 - x * x);
                    Pnm[3] = 15.0 * s * s * s;
                    break;
                case 4 :
                    Pnm[0] = (35.0 * x * x * x * x - 30.0 * x * x + 3.0) / 8.0;
                    Pnm[1] = 5.0 * (7.0 * x * x - 3.0) * x * s / 2.0;
                    Pnm[2] = 15.0 * (7.0 * x * x - 1.0) * (1.0 - x * x) / 2.0;
                    Pnm[3] = 105.0 * s * s * s * x;
                    Pnm[4] = 105.0 * s * s * s * s;
                    break;
                case 5 :
                    Pnm[0] = x * (63.0 * x * x * x * x - 70.0 * x * x + 15.0) / 8.0;
                    Pnm[1] = 15.0 * s * (21.0 * x * x * x * x - 14.0 * x * x + 1.0) / 8.0;
                    Pnm[2] = 105.0 * x * (3.0 * x * x - 1.0) * (1.0 - x * x) / 2.0;
                    Pnm[3] = 105.0 * s * s * s * (9.0 * x * x - 1.0) / 2.0;
                    Pnm[4] = 945.0 * s * s * s * s * x;
                    Pnm[5] = 945.0 * s * s * s * s * s;
                    break;
                case 6 :
                    Pnm[0] = (x * x * (x * x * (231.0 * x * x - 315.0) + 105.0) - 5.0) / 16.0;
                    Pnm[1] = 21.0 * x * (x * x * (33.0 * x * x - 30.0) + 5.0) * s / 8.0;
                    Pnm[2] = 105.0 * s * s * (x * x * (33.0 * x * x - 18.0) + 1.0) / 8.0;
                    Pnm[3] = 315.0 * (11.0 * x * x - 3.0) * x * s * s * s / 2.0;
                    Pnm[4] = 945.0 * s * s * s * s * (11.0 * x * x - 1.0) / 2.0;
                    Pnm[5] = 10395.0*x*s*s*s*s*s;
                    Pnm[6] = 10395.0*s*s*s*s*s*s;
                    break;
                case 7 :
                    Pnm[0] = x * (x * x * (429.0 *Math.Pow(x,4) - 693.0 * x * x + 315.0) - 35.0) / 16.0;
                    Pnm[1] = 7.0 * s * (x * x * (429.0 * Math.Pow(x, 4) - 495.0 * x * x + 135.0) - 5.0) / 16.0;
                    Pnm[2] = 63.0 * x * s * s * (x * x * (143.0 * x * x - 110.0) + 15.0) / 8.0;
                    Pnm[3] = 315.0 * s * s * s * (x * x * (143.0 * x * x - 66.0) + 3.0) / 8.0;
                    Pnm[4] = 3465.0 * x * s * s * s * s * (13.0 * x * x - 3.0) / 2.0;
                    Pnm[5] = 10395.0 * Math.Pow(s,5) * (13.0 * x * x - 1.0) / 2.0;
                    Pnm[6] = 135135.0 * x * Math.Pow(s ,6);
                    Pnm[7] = 135135.0 * Math.Pow(s,7);
                    break;
                default :
                    ALFS[0] = (x * x * (x * x * (x * x * (6435.0 * x * x - 12012.0) + 6930.0) - 1260.0) + 35) / 128.0;
                    ALFS[1] = 9.0 * x * s * (x * x * (x * x * (715.0 * x * x - 1001.0) + 385.0) - 35.0) / 16.0;
                    ALFS[2] = 315.0 * s * s * (x * x * (x * x * (143.0 * x * x - 143.0) + 33.0) - 1.0) / 16.0;
                    ALFS[3] = 3465.0 * x * s * s * s * (x * x * (39.0 * x * x - 26.0) + 3.0) / 8.0;
                    ALFS[4] = 10395.0 * s * s * s * s * (65.0 * Math.Pow(x, 4) - 26.0 * x * x + 1.0) / 8.0;
                    ALFS[5] = 135135.0 * x * Math.Pow(s ,5) * (5.0 * x * x - 1.0) / 2.0;
                    ALFS[6] = 135135.0 * Math.Pow(s ,6) * (15.0 * x * x - 1.0) / 2.0;
                    ALFS[7] = 2027025.0 * x * s * s * s * s * s * s * s;
                    ALFS[8] = 2027025.0 * s * s * s * s * s * s * s * s;

                    if (n == 8)
                    {
                        Pnm = ALFS;
                    }
                    else
                    {
                        BLFS[0] = x * (x * x * (429.0 *Math.Pow(x,4) - 693.0 * x * x + 315.0) - 35.0) / 16.0;
                        BLFS[1] = 7.0 * s * (x * x * (429.0 * Math.Pow(x, 4) - 495.0 * x * x + 135.0) - 5.0) / 16.0;
                        BLFS[2] = 63.0 * x * s * s * (x * x * (143.0 * x * x - 110.0) + 15.0) / 8.0;
                        BLFS[3] = 315.0 * s * s * s * (x * x * (143.0 * x * x - 66.0) + 3.0) / 8.0;
                        BLFS[4] = 3465.0 * x * s * s * s * s * (13.0 * x * x - 3.0) / 2.0;
                        BLFS[5] = 10395.0 * Math.Pow(s,5) * (13.0 * x * x - 1.0) / 2.0;
                        BLFS[6] = 135135.0 * x * Math.Pow(s, 6);
                        BLFS[7] = 135135.0 * Math.Pow(s, 7);

                        int l = 8;
                        while (l <n)
                        {
                            Pnm[0] = Legendre.LegendreF(l+1 , x);
                            int k = 1;
                            while (k <= l + 1)
                            {
                               if (k <= l - 1)
                                {
                                    Pnm[k] = ((2 * l + 1) * x * ALFS[k] - (l + k) * BLFS[k]) / (l - k + 1);
                                }
                                               
                               else if(k > l - 1)
                                {
                                    Pnm[k] = -((l - k + 2) * x * Pnm[k - 1] - (l + k) * ALFS[k - 1]) / s;
                                }
                                k = k + 1; 
                            }
                            Array.Copy(ALFS, BLFS, ALFS.Length);
                            Array.Copy(Pnm, ALFS, Pnm.Length);
                            l = l + 1;
                        }
                       
                     }
                break; 
            }   
          return Pnm;   
        }


        public static double[] ALegendFnmTest(int n, double x)
        {
            double[] Pnm = new double[n + 1];
            double[] ALFS = new double[n + 1];
            double[] BLFS = new double[n + 1];

            double s = Math.Sqrt(1.0 - x * x);

            switch (n)
            {
                case 0:                    
                    Pnm[0] = 1.0;
                    break;
                case 1:                   
                    Pnm[0] = x;
                    Pnm[1] = s;
                    break;
               default:                   
                    ALFS[0] = (3.0 * x * x - 1.0) / 2.0;
                    ALFS[1] = 3.0 * s * x;
                    ALFS[2] = 3.0 * (1.0 - x * x);

                    if (n == 2)
                    {
                        Pnm = ALFS;
                    }
                    else 
                    {                       
                        BLFS[0] = x;
                        BLFS[1] = s;

                        int l = 2;
                        while (l < n)
                        {
                            Pnm[0] = Legendre.LegendreF(l + 1, x);
                            int k = 1;
                            while (k <= l + 1)
                            {
                                if (k <= l - 1)
                                {                                    
                                    Pnm[k] = ((2 * l + 1) * x * ALFS[k] - (l + k) * BLFS[k]) / (l - k + 1);                                    
                                }

                                else if (k > l - 1)
                                {                                   
                                    Pnm[k] = -((l - k + 2) * x * Pnm[k - 1] - (l + k) * ALFS[k - 1]) / s;
                                }
                                k = k + 1;
                            }                            
                            Array.Copy(ALFS, BLFS, ALFS.Length);                           
                            Array.Copy(Pnm, ALFS, Pnm.Length);                         
                            l = l + 1;
                        }
                    }
               break; 
            }             
            return Pnm;
        }


        public static double[] FulNormALegenSerie(int n, double x)
        {
                       
            double[] Pnm = new double[n + 1];
            double[] BLFS = new double[n + 1];
            double[] ALFS = new double[n + 1];

            double s = Math.Sqrt(1.0 - x * x);

            if (n == 0)
            {                
                Pnm[0] = 1.0;
                return Pnm;
            }
            if (n == 1)
            {
                Pnm[0] = x * Math.Sqrt(3.0);
                Pnm[1] = s * Math.Sqrt(3.0);                
                return Pnm;
            }            
            else if(n == 2)
            {
                Pnm[0] = Math.Sqrt(5.0) * LegendreF(2, x);
                Pnm[1] = 3.0 * s * x * Math.Sqrt(10.0 / 6.0);
                Pnm[2] = 3.0 * (1.0 - x * x) * Math.Sqrt(10.0 / 24.0);                
                return Pnm;
            }              
            else
            {
                   BLFS[0] = x * Math.Sqrt(3.0);
                   BLFS[1] = s * Math.Sqrt(3.0);                

                  ALFS[0] = Math.Sqrt(5.0) * LegendreF(2, x);
                  ALFS[1] = 3.0 * s * x * Math.Sqrt(10.0 / 6.0);
                  ALFS[2] = 3.0 * (1.0 - x * x) * Math.Sqrt(10.0 / 24.0);                 
                              
                Pnm[0] = Math.Sqrt(2.0 * n + 1.0) * LegendreF(n, x);
                int l = 3;
                while (l <= n)
                {                    
                    int k = 1; 
                    while (k <= l)
                    {
                        if (k < l - 1)
                        {
                            double Anm = Math.Sqrt((2.0 * l - 1.0) * (2.0 * l + 1.0)/ 
                                ((l - k) * (l + k)));
                            double Bnm = -Math.Sqrt(((l - 1.0) * (l - 1.0) - k * k) 
                                / (4.0 * (l - 1.0) * (l - 1.0) - 1.0));
                            Pnm[k] = Anm * (x * ALFS[k] + Bnm * BLFS[k]);

                        }
                        else  if(k >= l - 1)
                        {
                            if (k == l - 1)
                            {
                                Pnm[k]= Math.Sqrt(2.0 * k + 3.0) * x * ALFS[l - 1];
                            }                                    
                            else if(k == l)
                            {
                                Pnm[k] = s * Math.Sqrt(1.0 + 1.0 / (2.0 * k)) * ALFS[l - 1];
                            }  
                        }
                        k = k + 1;                  
                    }
                    Array.Copy(ALFS, BLFS, ALFS.Length);
                    Array.Copy(Pnm, ALFS, Pnm.Length);
                    l = l + 1;
                }
                return Pnm;
            } 
        }

        public static double[] ellipsoidal (double R , double latitute)
        {
            double Rlat = latitute * DEGRAD;
            double x = Math.Sin(Rlat);
            double x2 = x * x;
            double w = Math.Sqrt(1.0 - e2 * x2);
            double rp = Math.Sqrt(1.0 + e2 * (e2 - 2.0) * x2) * R / w;
            double Rn = R / rp;
            double rep = R/Math.Sqrt(1.0 + ep2 * x2) ;
            double Ren = R / rep;

            double[] ellipsoid =  {rp , Rn , rep , Ren };

            return ellipsoid;
        }

        public static double gravitation(double latitude, double longitude, int method , string fileName)
        {

          double Rlam = longitude * DEGRAD;
          double x = 0.0; 

           double Sum1 = 0.0;
           double Tsum = 0.0;
           double G = 0.0;
           double R = 0.0;
           double r = 0.0;
            double Rn = 0.0;
           int n = 0;

            string line = "";
            try
            {
                using (System.IO.StreamReader sr = new System.IO.StreamReader(fileName)) 
                {
                    while ((line = sr.ReadLine()) != null)
                    {
                        if (line.StartsWith("earth_gravity_constant"))
                        {
                            string[] ar = line.Split(' ', 2);
                            ar[1] = ar[1].Trim();
                            if (ar[1].Contains('D')) 
                                {
                                   ar[1] = ar[1].Replace("D","E"); 
                                }
                            G = Convert.ToDouble(ar[1]);
                            Console.WriteLine(" earth_gravity_constant= {0}", G);
                        }
                        else if (line.StartsWith("radius"))
                        {
                            string[] ar = line.Split(' ', 2);
                            ar[1] = ar[1].Trim();
                            if (ar[1].Contains('D'))
                            {
                                ar[1] = ar[1].Replace("D", "E");
                            }
                            R = Convert.ToDouble(ar[1]);                            
                            Console.WriteLine(" radius= {0}", R);
                        }
                        else if (line.StartsWith("max_degree"))
                        {
                            string[] ar = line.Split(' ', 2);
                            n = Convert.ToInt32(ar[1].Trim());
                            Console.WriteLine(" max_degree= {0}", n);
                            break;
                        }
                    }

                    double[] ellipsoid = ellipsoidal(R, latitude);

                    if (method == 0)
                    {
                        r = ellipsoid[0];
                        Rn = ellipsoid[1];
                    }
                    else if (method == 1)
                    {
                        r = ellipsoid[2];
                        Rn = ellipsoid[3];                    
                    }
                    
                    x = Math.Sin(latitude * DEGRAD); 
                    double[] Pn = new double[n+1];

                    while ((line = sr.ReadLine()) != null)
                    { 
                        if (line.StartsWith("gfc"))
                        {
                            string ar = line.Split(' ',2)[1];                            
                            int L  = Convert.ToInt32(ar.Trim().Split(' ', 2)[0]);
                            string br = ar.Trim().Split(' ', 2)[1];
                            int M = Convert.ToInt32(br.Trim().Split(' ',2)[0]);
                            string cr = br.Trim().Split(' ', 2)[1];
                            double C = Convert.ToDouble(cr.Trim().Split(' ', 2)[0]);
                            string dr = cr.Trim().Split(' ', 2)[1];
                            double S = Convert.ToDouble(dr.Trim().Split(' ', 2)[0]); 
                            
                            if (M == 0)
                            {                                
                                Pn = FulNormALegenSerie(L, x);
                                Sum1 = Pn[M] * C ;                               
                            }
                            else 
                               {                                                              
                                Sum1 = Sum1 + Pn[M] * (C * Math.Cos(M * Rlam) + 
                                    S * Math.Sin(M * Rlam)) ;
                               }
                             if (M == L)
                               {                                
                                Tsum = Tsum + Sum1 * (L-1)*Math.Pow(Rn,L);                                
                               }                            
                        }                        
                    }
                }                
            }
            catch (Exception e)
            {
                Console.WriteLine("Error in Reading File:");
                Console.WriteLine(e.Message);
            }
            double Coef = G  / (r * r);
            return Tsum * Coef;
        }

        public static double normalGravity(double latitute)
        {
            double Rlat = latitute * DEGRAD;
            double x = Math.Sin(Rlat);

            double Term = 0.0;
            double e4 = e2 * e2;
            double e6 = e4 * e2;

            double[] a2n = { 0.50 * e2 + k ,
                            3.0 * e4 / 8.0 + 0.50 * e2 * k ,
                            5.0 * e6 / 16.0 + 3.0 * e4 * k / 8.0 ,
                            35.0 * e4 * e4 / 128.0 + 5.0 * e6 * k / 16.0 };

            for(int l = 1; l<=4; l++)
            {
                Term = Term + a2n[l - 1] * Math.Pow(x, 2 * l);
            }

                return  9.78032677150 * (1.0 + Term);

        }

        public static double adaptedGravity(double latitute, double height)
        {

            double ga = 9.78032677150;

            double Rlat = latitute * DEGRAD;
            double x = Math.Sin(Rlat);            
            double x2 = x * x;

            double adaptedG = ga * (1 + x2 * (0.00527904140 + x2 * (0.00002327180 +
              x2 * (0.00000012620 + 0.00000000070 * x2)))) - height * (0.03087798E-4 -
                    x2 * (0.0000439E-4 + 0.00000020E-4 * x2)) - height * height * (-0.00007265E-8
                     + 0.00000021E-8 * x2);

            return adaptedG;
        }

        public static double Calcgravity(double latitute, int method)
        {

            double[] C_Coef = new double[5];
            double[] ellipsoid = ellipsoidal(a, latitute);
            double r = 0.0;
            double Rn = 0.0;

            if (method == 0)
            {
                r = ellipsoid[0];
                Rn = ellipsoid[1];
                C_Coef[0] = 1.0000055695713940820;
                C_Coef[1] = -0.0010812989696786500;
                C_Coef[2] = 0.0000004713225759810;
                C_Coef[3] = -0.0000000009225623140;
                C_Coef[4] = 0.0000000000104590330;
            }
            else if (method == 1)
            {
                r = ellipsoid[2];
                Rn = ellipsoid[3];
                C_Coef[0] = 0.9999935514460012160;
                C_Coef[1] = -0.0010841619139521870;
                C_Coef[2] = 0.0000045717325681810;
                C_Coef[3] = -0.0000000170436326690;
                C_Coef[4] = 0.0000000000817350280;
            }
                                  
            double Rphi = latitute * DEGRAD;
            double x = Math.Sin(Rphi);
            double s = Math.Cos(Rphi);
            double x2 = x * x;
            double s2 = s * s;
            double Omega2 = Omega * Omega;

            double Sigma = 0.0;

             for(int j = 0; j<= 4; j++)
            {
                Sigma = Sigma + GM * ((2 * j + 1) * C_Coef[j] * LegendreF(2 * j, x) * Math.Pow(Rn, 2 * j + 1));
            }

            return Sigma / (r * a) - r * Omega2 * s2;
        }
    }  
}
