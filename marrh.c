#include <stdio.h>                  /*  Marr-Hildreth.c  (or marrh.c) */
#include <math.h>
#include <stdlib.h>
#define  PICSIZE 256
#define  MAXMASK 100

int    pic[PICSIZE][PICSIZE];
double outpic1[PICSIZE][PICSIZE];
double outpic2[PICSIZE][PICSIZE];
int    edgeflag[PICSIZE][PICSIZE];
double maskx[MAXMASK][MAXMASK];
double masky[MAXMASK][MAXMASK];
double convx[PICSIZE][PICSIZE];
double convy[PICSIZE][PICSIZE];
double ival[PICSIZE][PICSIZE];
double peaks[PICSIZE][PICSIZE];

int main(argc,argv)
int argc;
char **argv;
{
    // Local vars
    int     i,j,p,q,s,t,mr,centx,centy;
    double  maskvalx,maskvaly,sumx,sumy,sig,maxival,minival,maxval,ZEROTOL, slope;
    FILE    *fo1, *fo2, *fp1, *fopen();
    char    *foobar;
    
    // Input file
    fp1=fopen("garb34.pgm","rb");
    
    // Output files
    fo1=fopen("out1.pgm","wb");
    fo2=fopen("out2.pgm","wb");
    
    // Sigma value
    sig = 1.0;
    
    // Zero tolerance
    ZEROTOL = atof(foobar);
    
    mr = (int)(sig * 3);
    centx = (MAXMASK / 2);
    centy = (MAXMASK / 2);
    
    // Remove file header on input
    char tmp[20];
    for (i=1;i<=4;i++) {
        fscanf(fp1,"%s", tmp);
    }
    
    // Get values from input picture
    for (i=0;i<256;i++) {
        for (j=0;j<256;j++) {
            pic[i][j]  =  getc(fp1);
            pic[i][j]  &= 0377;
        }
    }
    
    // Get derivitaves
    for (p=-mr;p<=mr;p++) {
        for (q=-mr;q<=mr;q++) {
            maskvalx = p * (exp(-1 * ( ((p*p)+(q*q)) / (2 * (sig*sig) ))));
            maskvaly = q * (exp(-1 * ( ((p*p)+(q*q)) / (2 * (sig*sig) ))));
            (maskx[p+centy][q+centx]) = maskvalx;
            (masky[p+centy][q+centx]) = maskvaly;
        }
    }
    
    // Convolution code
    for (i=mr;i<=255;i++) {
        for (j=mr;j<=255;j++) {
            sumx = 0;
            sumy = 0;
            for (p=-mr;p<=mr;p++) {
                for (q=-mr;q<=mr;q++) {
                    sumx += pic[i+p][j+q] * maskx[p+centy][q+centx];
                    sumy += pic[i+p][j+q] * masky[p+centy][q+centx];
                }
            }
            outpic1[i][j] = sumx;
            outpic2[i][j] = sumy;
            convx[i][j] = sumx;
            convy[i][j] = sumy;
        }
    }
    
    // Format headings of files
    fprintf(fo1, "P5\n");
    fprintf(fo1, "%d %d\n", 256, 256);
    fprintf(fo1, "255\n");
    
    
    // Get ival and find maxival
    maxival = 0;
    for (i=mr; i<256-mr; i++) {
        for (j=mr; j<256-mr; j++) {
            ival[i][j]=sqrt((double)((outpic1[i][j] * outpic1[i][j]) +
                                     (outpic2[i][j] * outpic2[i][j])));
            if (ival[i][j] > maxival)
                maxival = ival[i][j];
        }
    }
    
    // Magnitude algorithm
    for (i=0; i<256; i++) {
        for (j=0; j<256; j++) {
            ival[i][j] = (ival[i][j] / maxival) * 255;
            fprintf(fo1,"%c",(char)((int)(ival[i][j])));
        }
    }
    

    //Thresholds
    for (i=0;i<256-mr;i++) {
        for (j=0;j<256-mr;j++) {
            //Low
            if(ival[i][j] > maxival)
                    fprintf(fo1,"%c",(char)((int)(255)));
            else
                    fprintf(fo1,"%c",(char)((int)(0)));
        }
    }
    
    //------------------------------ END OF PART ONE --------------------------

    // Format headings of file 2
    fprintf(fo2, "P5\n");
    fprintf(fo2, "%d %d\n", 256, 256);
    fprintf(fo2, "255\n");


    // Loop
    for (i=mr;i<256-mr;i++) {
      for (j=mr;j<256-mr;j++) {

        // Make sure we dont divide by 0
        if (ival[i][j] == 0.0)
            ival[i][j] = 0.00001;

        // Get slope
        slope = convy[i][j] / convx[i][j];
        
        // Initialize
        peaks[i][j] = 0;

        // Which cone?
        if ( (slope <= 0.4142) && (slope > -0.4142) ) {
            if (ival[i][j] > ival[i][j-1] && ival[i][j] > ival[i][j+1])
                peaks[i][j] = 255;
        }
        else if ( (slope <= 2.4142) && (slope > 0.4142) ) {
            if (ival[i][j] > ival[i-1][j-1] && ival[i][j] > ival[i+1][j+1])
                peaks[i][j] = 255;
        }
        else if ( (slope <= -0.4142) && (slope > -2.4142) ) {
            if (ival[i][j] > ival[i+1][j-1] && ival[i][j] > ival[i-1][j+1])
                peaks[i][j] = 255;
        }
        else {
            if (ival[i][j] > ival[i-1][j] && ival[i][j] > ival[i+1][j])
                peaks[i][j] = 255;
        }
        }
    }

    //Print image
    for (i=0;i<256;i++) { 
        for (j=0;j<256;j++) {  
                if(peaks[i][j] == 255)
                    fprintf(fo2,"%c",(char)((int)(255))); //(peaks[i][j])));
                else
                    fprintf(fo2,"%c",(char)((int)(0)));
                //fprintf(fo2,"%c",(char)((int)(peaks[i][j])));
        }
    }

    


    // ----------------------------- END OF PART TWO --------------------------
    //  // Max vals
    //  maxval  = 0;
    //  maxival = 0;
    //  minival = 255;
    
    // // Find min/max vals
    // for (i=mr;i<256-mr;i++) {
    //   for (j=mr;j<256-mr;j++) {
    //       if (outpic1[i][j] == 0.0)
    //           outpic1[i][j] = 0.00001;
    //       if (outpic1[i][j] > maxival )
    //           maxival = outpic1[i][j];
    //       if (outpic1[i][j] < minival)
    //           minival = outpic1[i][j];
    //   }
    // }
    
    // if (fabs(maxival) > fabs(minival))
    //     maxval = fabs(maxival);
    // else
    //     maxval = fabs(minival);

    // for (i=0;i<256;i++) { 
    //     for (j=0;j<256;j++) {
    //         outpic1[i][j] = ((((outpic1[i][j]) / maxval) + 1) * 127);
    //         fprintf(fo1,"%c",(char)((int)(outpic1[i][j])));
    //     }
    // }
    
    // for (i=mr;i<=255-mr;i++) {
    //     for (j=mr;j<=255-mr;j++) {

    //         outpic2[i][j] = 0;
    //         if (conv[i][j] > ZEROTOL)
    //         {
    //         for (p=-1;p<=1;p++)
    //         {
    //             for (q=-1;q<=1;q++)
    //             {
    //             if (conv[i+p][j+q] < -ZEROTOL)
    //             {
    //                 outpic2[i][j] = 255;
    //             }
    //             }
    //         }
    //         }
    //         else if ((fabs)(conv[i][j]) < ZEROTOL)
    //         {
    //                 if (((conv[i+1][j] > ZEROTOL) &&
    //                     (conv[i-1][j] < -ZEROTOL))   ||
    //                     ((conv[i+1][j] < -ZEROTOL) &&
    //                     (conv[i-1][j] > ZEROTOL)))
    //                 {
    //                 outpic2[i][j] = 255;
    //                 }
    //                 else if (((conv[i][j+1] > ZEROTOL) &&
    //                         (conv[i][j-1] < -ZEROTOL))   ||
    //                         ((conv[i][j+1] < -ZEROTOL) &&
    //                         (conv[i][j-1] > ZEROTOL)))
    //                 {
    //                 outpic2[i][j] = 255;
    //                 }
    //                 else if (((conv[i+1][j+1] > ZEROTOL) &&
    //                         (conv[i-1][j-1] < -ZEROTOL))   ||
    //                         ((conv[i+1][j+1] < -ZEROTOL) &&
    //                         (conv[i-1][j-1] > ZEROTOL)))
    //                 {
    //                 outpic2[i][j] = 255;
    //                 }
    //                 else if (((conv[i+1][j-1] > ZEROTOL) &&
    //                         (conv[i-1][j+1] < -ZEROTOL))   ||
    //                         ((conv[i+1][j-1] < -ZEROTOL) &&
    //                         (conv[i-1][j+1] > ZEROTOL)))
    //                 {
    //                 outpic2[i][j] = 255;
    //                 }
    //         }
    //     }
    // }
    
        // for (i=0;i<256;i++)
        // { for (j=0;j<256;j++)
        //   {  if (outpic2[i][j] == 255) outpic2[i][j]=0;
        //       else outpic2[i][j]=255;
        //       fprintf(fo2,"%c",(char)((int)(outpic2[i][j])));
        //   }
        // }
 }
