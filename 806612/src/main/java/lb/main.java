/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lb;

import java.awt.* ;
import javax.swing.* ;
import java.util.Scanner;

import com.aparapi.*;
import com.aparapi.device.Device ;

public class main {

    final static int maxIter = 200000;

    final static float Re = 220.0F;  // Reynolds number

    final static int nx = 520, ny = 180;  // Lattice dimensions
    final static int q = 9;  // population

    final static float uLB = 0.06F;  // Velocity in lattice units

    // useful constants for weights

    final static float W0 = 4.0F / 9;
    final static float W1 = 1.0F / 9;
    final static float W2 = 1.0F / 36;

    final static int CELL_SIZE = 2;
    final static int OUTPUT_FREQ = 100;

    static int [] [] c = new int [q] [2];  // Lattice velocities
    static double [] t = new double [q];  // Lattice weights
    
    final static float [] u = new float [nx * ny * 2];

    public static void main(String args []) throws Exception {

        float [] u = new float [nx * ny * 2];
        
        float ly = ny - 1.0F;

        float cx = nx/4.0F, cy = ny/2.0F, r = ny/9.0F;
        // Coordinates, size of cylinder.

        float nulb = uLB * r / Re;
        float omega = 1.0F / (3 * nulb + 0.5F);  // Relaxation parameter

        long startTime = System.currentTimeMillis();

        long macroTime = 0;
        long collisionTime = 0;
        long streamingTime = 0;

        int [] vStates = new int [] {0, -1, +1};
        int pos = 0;
        for(int i : vStates)
            for(int j : vStates) {
                int [] cEl = c [pos++];
                cEl [0] = i;
                cEl [1] = j;
            }

        t [0] = W0;
        for(int i = 1; i < q; i++) {
            int [] cEl = c [i];
            if(cEl [0] == 0 || cEl [1] == 0) {
                t [i] = W1;
            } else {
                t [i] = W2;
            }
        }

        int [] noslip = new int [q];  // index in c of negative velocity state
        for(int i = 0; i < q; i++) {
            int [] cEl = c [i];
            for(int j = 0; j < q; j++) {
                int [] cElj = c [j];
                if(cElj [0] == -cEl [0] && cElj [1] == -cEl [1]) {
                    noslip [i] = j;
                }
            }
        }

        int [] i1 = new int [3], i2 = new int [3], i3 = new int [3];
        int i1pos = 0, i2pos = 0, i3pos = 0;
        for(int i = 0; i < q; i++) {
            int cElX = c [i] [0];
            if (cElX < 0) {
                i1 [i1pos++] = i;
            } else if (cElX == 0) {
                i2 [i2pos++] = i;
            } else {
                i3 [i3pos++] = i;
            }
        }

        // Cylindrical (0) or Airfoil (1) obstacle
        Scanner sc = new Scanner(System.in);
        System.out.println("Objects available are:\n0 - Cylinder\n1 - Airfoil\n"
                                + "Enter an object from above as an integer: ");
        int inputObject;
            while (true){
                if (sc.hasNextInt()){
                    inputObject = sc.nextInt(); // Assign the next integer to a variable
                    if (inputObject  <= 255  && inputObject >= 0){ // Check if integer meets condition
                        break; // Condition met, break out of loop
                    } else if (inputObject > 255) {
                        System.out.println("Number cannot be over 255\n"
                                + "Please enter an integer between 0 and 255: ");
                    } else {
                        System.out.println("Number cannot be negative!\n"
                                + "Please enter a number between 0 and 255 ");
                    }
                } else {
                    System.out.println("Value entered is not an integer!\n"
                                + "Please enter an integer between 0 and 255: ");
                    sc.next();
                }
            }
        int object = inputObject;
        boolean [] obstacle = new boolean [nx * ny];
        double r2 = r * r;
        if (object == 1) {
            System.out.println("Object is 1 - Airfoil");
            Scanner airfoilNum = new Scanner(System.in);
            System.out.println("A 4-digit NACA number XYZZ is needed "
                                + "to calculate shape of airfoil.\n"
                                + "Enter a four digit integer not divisible by 100: ");
            int inputNaca;
            while (true){
                if (airfoilNum.hasNextInt()){
                    inputNaca = airfoilNum.nextInt(); // Assign the next integer to a variable
                    if (inputNaca  <= 9999  && inputNaca >= 1){ // Check if integer meets condition
                        if (inputNaca % 100 == 0) {
                            System.out.println("Last two digits cannot be 00!\n"
                                + "Please enter an integer that is not divisible by 100: ");
                        } else {
                            break; // Condition met, break out of loop
                        }
                    } else if (inputNaca > 9999) {
                        System.out.println("Number cannot be over 4 digits!\n"
                                + "Please enter an integer with at most 4 digits: ");
                    } else if (inputNaca == 0) {
                        System.out.println("Number cannot be 0!\n"
                                + "Please enter a positive 4-digit number: ");
                    } else {
                        System.out.println("Number cannot be negative!\n"
                                + "Please enter a positive 4-digit number: ");
                    }
                } else {
                    System.out.println("Value entered is not an integer!\n"
                                + "Please enter a 4-digit integer: ");
                    sc.next();
                }
            }
            int nacaNum = inputNaca;
            System.out.println("NACA number entered is " + nacaNum);
            double nacax = (int) (nacaNum / 1000);
            double nacay = (int) ((nacaNum - (nacax * 1000)) / 100);
            double nacazz = (int) (nacaNum - (nacax * 1000) - (nacay * 100));
            System.out.println("X is: " + (int) nacax + ", Y is: " + (int) nacay
                                + ", and ZZ is: " + (int) nacazz);
            double [] ys = new double [2];
            for(int k = 0; k <= (q*20); k++) {
                double x = (double) k / (q*20);
                nacaxyzz(ys, x, nacax, nacay, nacazz);
                int i = (nx - (q*20)) / 4  + k;
                for(int j = 0; j < ny; j++) {
                    obstacle [ny * i + j] = (double) j > (ny/2 - (q*20) * ys [1]) &&
                                        (double) j < (ny/2 - (q*20) * ys [0]);
                }
            }
        } else {  
            if (object == 0) {
                System.out.println("Object is 0 - Cylinder"); 
            } else {
                System.out.println("The number entered has no assigned object!\n"
                                            + "Object defaults to 0 - Cylinder");
            }
            
            for(int i = 0; i < nx; i++) {
                for(int j = 0; j < ny; j++) {
                    double dx = i - cx;
                    double dy = j - cy;
                    obstacle [ny * i + j] = (dx * dx + dy * dy) < r2;
                }
            }
            
        }

        // Velocity inlet with perturbation
        float [] vel = new float [ny * 2];
        for(int j = 0; j < ny; j++) {
            vel [j * 2 + 0] = (float) (uLB * (1 + 1E-4 * Math.sin((j/ly) * 2 * Math.PI)));
        }

        float [] fin = new float [nx * ny * q];
        for(int i = 0; i < nx; i++) {
            for(int j = 0; j < ny; j++) {
                if(i == 0) {
                    equilibrium(fin, (i * ny + j) * q, 1.0F, vel [j * 2 + 0], vel [j * 2 + 1]);
                }
                else {
                    equilibrium(fin, (i * ny + j) * q, 1.0F, 0.0F, 0.0F);
                }
            }
        }

        float [] fout = new float [nx * ny * q];
        float [] rho = new float [nx * ny];
        
        Display display = new Display(u);
        
        final int[] time = new int[]{0};
        Device device = Device.best();
        Range range = device.createRange2D(nx, ny);
        
        float [] feq = new float [nx * ny * q];
        
        Kernel kernelMacroCollision = new Kernel() {
            @Override
            public void run() {
                int i = getGlobalId(0);
                int j = getGlobalId(1);
                int base_ij = (i * ny + j) * q;
                int base_uij = (i * ny + j) * 2;
                if(i > 0) {
                    float sum = fin [base_ij + 0] + fin [base_ij + 1] + fin [base_ij + 2] +
                                fin [base_ij + 3] + fin [base_ij + 4] + fin [base_ij + 5] +
                                fin [base_ij + 6] + fin [base_ij + 7] + fin [base_ij + 8];

                    float sum0 = - fin [base_ij + 3] - fin [base_ij + 4] - fin [base_ij + 5]
                                 + fin [base_ij + 6] + fin [base_ij + 7] + fin [base_ij + 8];

                    float sum1 = - fin [base_ij + 1] + fin [base_ij + 2] - fin [base_ij + 4]
                                 + fin [base_ij + 5] - fin [base_ij + 7] + fin [base_ij + 8];

                    rho [ny * i + j] = sum;
                    if(sum > 0) {
                        u [base_uij + 0] = sum0 / sum;
                        u [base_uij + 1] = sum1 / sum;
                    }
                } else {
                    u [base_uij + 0] = (float) (vel [j * 2 + 0]);
                    u [base_uij + 1] = (float) (vel [j * 2 + 1]);
                          
                    float sum2 = 0;
                    for(int k = 0; k < 3; k++) {
                        int d = i2 [k];
                        sum2 += fin [base_ij + d];
                    }
                    float sum1 = 0;
                    for(int k = 0; k < 3; k++) {
                        int d = i1 [k];
                        sum1 += fin [base_ij + d];
                    }
                    rho [0 + j] = (float) (1/(1 - u [base_uij + 0]) * (sum2 + 2 * sum1));
                }
                // Collision step.
                if(obstacle [ny * i + j]) {
                    for(int d = 0; d < q; d++) {
                        fout [base_ij + d] = fin [base_ij + noslip [d]];
                    }
                } else {
                    equilibrium(feq, base_ij, rho [ny * i + j], u [((i) * ny + j) * 2 + 0], u [((i) * ny + j) * 2 + 1]);
                    if(i == 0) {
                        for(int p = 0; p < 3; p++) {
                            fin [base_ij + i3 [p]] = feq [base_ij + i3 [p]];
                        }
                    }
                    fout [base_ij + 0] = (fin [base_ij + 0] -
                                omega * (fin [base_ij + 0] - feq [base_ij + 0]));
                    fout [base_ij + 1] = (fin [base_ij + 1] -
                                omega * (fin [base_ij + 1] - feq [base_ij + 1]));
                    fout [base_ij + 2] = (fin [base_ij + 2] -
                                omega * (fin [base_ij + 2] - feq [base_ij + 2]));
                    fout [base_ij + 3] = (fin [base_ij + 3] -
                                omega * (fin [base_ij + 3] - feq [base_ij + 3]));
                    fout [base_ij + 4] = (fin [base_ij + 4] -
                                omega * (fin [base_ij + 4] - feq [base_ij + 4]));
                    fout [base_ij + 5] = (fin [base_ij + 5] -
                                omega * (fin [base_ij + 5] - feq [base_ij + 5]));
                    fout [base_ij + 6] = (fin [base_ij + 6] -
                                omega * (fin [base_ij + 6] - feq [base_ij + 6]));
                    fout [base_ij + 7] = (fin [base_ij + 7] -
                                omega * (fin [base_ij + 7] - feq [base_ij + 7]));
                    fout [base_ij + 8] = (fin [base_ij + 8] -
                                omega * (fin [base_ij + 8] - feq [base_ij + 8]));
                }
            }
        };
            
        // Streaming step.
        Kernel kernelStream = new Kernel() {
            @Override
            public void run() {
                int i = getGlobalId(0);
                
                int iP1 = (i + 1) % nx;
                int iM1 = (i - 1 + nx) % nx;

                int base_i = i * ny * q;

                int base_iM1 = iM1 * ny * q;

                int base_iP1 = iP1 * ny * q;


                int j = getGlobalId(1);
                int base_ij = (i * ny + j) * q;

                int jP1 = (j + 1) % ny;
                int jM1 = (j - 1 + ny) % ny;
                    
                int offs_j = j * q ;
                int offs_jP1 = jP1 * q ;
                int offs_jM1 = jM1 * q ;  // etc, maybe

                fin [base_i + offs_j + 0] = fout [base_ij + 0];
                fin [base_i + offs_jM1 + 1] = fout [base_ij + 1];
                fin [base_i + offs_jP1 + 2] = fout [base_ij + 2];
                fin [base_iM1 + offs_j + 3] = fout [base_ij + 3];
                fin [base_iM1 + offs_jM1 + 4] = fout [base_ij + 4];
                fin [base_iM1 + offs_jP1 + 5] = fout [base_ij + 5];
                fin [base_iP1 + offs_j + 6] = fout [base_ij + 6];
                fin [base_iP1 + offs_jM1 + 7] = fout [base_ij + 7];
                fin [base_iP1 + offs_jP1 + 8] = fout [base_ij + 8];
            }
        };
         
        kernelMacroCollision.setExplicit(true);
        
        kernelMacroCollision.put(vel);
        
        kernelMacroCollision.put(feq);
        
        kernelStream.setExplicit(true);
        
        time[0]=0;
        
        for (int timeInt = 0; timeInt <= maxIter; timeInt++) {
            
            long time1 = System.currentTimeMillis();

            // Calculate macroscopic density and velocity
            kernelMacroCollision.put(fin);
            
            long time2 = System.currentTimeMillis();

            macroTime += (time2 - time1);

            //kernelMacroCollision.put(rho);
            //kernelMacroCollision.put(u);
            kernelMacroCollision.execute(range).get(fout);
            
            long time3 = System.currentTimeMillis();

            collisionTime += (time3 - time2);
            
            kernelStream.put(fout);
            kernelStream.execute(range).get(fin);
            
            long time4 = System.currentTimeMillis();

            streamingTime += (time4 - time3) ;
            
            // Right wall: outflow condition
            for(int d : i1) {
                for(int j = 0; j < ny; j++) {
                    fin [((nx - 1) * ny  + j) * q + d] = fin [((nx - 2) * ny + j) * q + d];
                }
            }
            
            if(time[0] % OUTPUT_FREQ == 0) { // change time[0] to timeInt for diagnosis purposes.
                //System.out.println("Time: " + timeInt);
                kernelMacroCollision.get(u);
                display.repaint();
            }
        }


        long endTime = System.currentTimeMillis();

        System.out.println("Calculation completed in " +
                            (endTime - startTime) + " milliseconds");

        System.out.println("Time to calculate macroscopic quantities: " +
                            macroTime + " milliseconds");
        System.out.println("Time for collision steps: " +
                            collisionTime + " milliseconds");
        System.out.println("Time for streaming steps: " +
                            streamingTime + " milliseconds");
        
        
        display.repaint();
    }


    static void equilibrium(float [] feq, int fbase, float rho,
                            float u0, float u1) {

        float usqr = u0 * u0 + u1 * u1;

        float u0Pu1 = u0 + u1;
        float u0Mu1 = u0 - u1;

        feq [fbase + 0] = (rho * W0 * (1.0F - 1.5F * usqr));
        feq [fbase + 1] = (rho * W1 * (1.0F - 3.0F * u1 + 4.5F *
                                                u1 * u1 - 1.5F * usqr));
        feq [fbase + 2] = (rho * W1 * (1.0F + 3.0F * u1 + 4.5F *
                                                u1 * u1 - 1.5F * usqr));
        feq [fbase + 3] = (rho * W1 * (1.0F - 3.0F * u0 + 4.5F *
                                                u0 * u0 - 1.5F * usqr));
        feq [fbase + 4] = (rho * W2 * (1.0F - 3.0F * u0Pu1 + 4.5F
                                                * u0Pu1 * u0Pu1 - 1.5F * usqr));
        feq [fbase + 5] = (rho * W2 * (1.0F - 3.0F * u0Mu1 + 4.5F
                                                * u0Mu1 * u0Mu1 - 1.5F * usqr));
        feq [fbase + 6] = (rho * W1 * (1.0F + 3.0F * u0 + 4.5F *
                                                u0 * u0 - 1.5F * usqr));
        feq [fbase + 7] = (rho * W2 * (1.0F + 3.0F * u0Mu1 +
                                            4.5F * u0Mu1 * u0Mu1 - 1.5F * usqr));
        feq [fbase + 8] = (rho * W2 * (1.0F + 3.0F * u0Pu1 +
                                            4.5F * u0Pu1 * u0Pu1 - 1.5F * usqr));
    }

    static class Display extends JPanel {
        float [] u ;

        Display(float [] u) {
            setPreferredSize(new Dimension(CELL_SIZE * nx, CELL_SIZE * ny));

            JFrame frame = new JFrame("LBM");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);

            this.u = u;            
        }

        public void paintComponent(Graphics g) {
            //double [] [] usqr = new double [nx] [ny];
            double [] usqr = new double [nx * ny];
            double uMax = Double.MIN_VALUE;
            double uMin = Double.MAX_VALUE;
            for(int i = 0; i < nx; i++) {
                for(int j = 0; j < ny; j++) {
                    double u0 = u [(i * ny + j) * 2 + 0];
                    double u1 = u [(i * ny + j) * 2 + 1];
                    double u2 = Math.sqrt(u0 * u0 + u1 * u1);
                    if(u2 < uMin) {
                        uMin = u2;
                    }
                    if(u2 > uMax) {
                        uMax = u2;
                    }
                    //usqr [i] [j] = u2;
                    usqr [ny * i + j] = u2;

                }
            }
            double diff = uMax - uMin;
            double norm = ((diff == 0.0) ? 0.0 : 1/diff);
            for(int i = 0; i < nx; i++) {
                for(int j = 0; j < ny; j++) {
                    //float f = (float) (norm * (usqr [i] [j] - uMin));
                    float f = (float) (norm * (usqr [ny * i + j] - uMin));
                    Color c = new Color(f, 0.0F, 1.0F - f);
                    g.setColor(c);
                    g.fillRect(CELL_SIZE * i, CELL_SIZE * j,
                                CELL_SIZE, CELL_SIZE);
                }
            }
            g.setColor(Color.WHITE);
            for(int i = 0; i < nx; i+=8) {
                for(int j = 0; j < ny; j+=8) {
                    int originX = CELL_SIZE * i + 5;
                    int originY = CELL_SIZE * j + 5;
                    g.drawOval(originX - 1, originY - 1, 3, 3);
                    g.drawLine(originX, originY,
                                originX + (int) (200 * u [(i * ny + j) * 2 + 0]),
                                originY + (int) (200 * u [(i * ny + j) * 2 + 1]));
                }
            }
        }
    }

    public final static double EPSILON = 0.000001;

    static void nacaxyzz(double [] ys, double x, double nacax, double nacay, double nacazz) {

        // https://en.wikipedia.org/wiki/NACA_airfoil


        // x argument here is actually x_U, x_L for upper, lower surface
        // respectively.

        // ys is two-element array containing y_L, y_U.
        
        for(int i = 0; i < 2; i++) {

            int sign = 2 * i - 1;  // -1, +1 yield y_L, y_U

            double xChord = x;
                    // xChord will be x/c in notation of Wiki page.

            // Solve for xChord in terms of x_U or x_L.
            double xOld;
            do {
                xOld = xChord;
                xChord = x + sign * y_t(xOld, nacazz) *
                                    Math.sin(theta(xOld, nacax, nacay));
            } while (Math.abs(xChord - xOld) > EPSILON);

            ys [i] = y_c(xChord, nacax, nacay) + sign * y_t(xChord, nacazz) *
                            Math.cos(theta(xChord, nacax, nacay));
        }
    }

    static double y_t(double x, double nacazz) {
        
        double t = nacazz / 100;  // (last two digits of 4415)
        double halfThickness = 5 * t * (0.2969 * Math.sqrt(x) + x *
                (-0.1260  + x * (-0.3516 + x * (0.2843 - x * 0.1015))));

        return halfThickness;
    }

    static double y_c(double x, double nacax, double nacay) {

        double m = nacax / 100;  // (first digit of 4415)
        double p = nacay / 10;   // (second digit of 4415)

        if(x < p) {
            return m * x * (2 * p - x) / (p * p);
        } else {
            return m * (1 - 2 * p + x * (2 * p - x)) / ((1 - p) * (1 - p));
        }
    }

    static double theta(double x, double nacax, double nacay) {

        double m = nacax / 100;  // (first digit of 4415)
        double p = nacay / 10;   // (second digit of 4415)

        double derivative;
        if(x < p) {
            derivative = 2 * m * (p - x) / (p * p);
        } else {
            derivative = 2 * m * (p - x) / ((1 - p) * (1 - p));
        }

        return Math.atan(derivative);
    }
}
