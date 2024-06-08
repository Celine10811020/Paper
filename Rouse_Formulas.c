#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct
{
    char Date[10];
    float DNALand;
    float DNA06m;
    float DNA15m;
    float Temp;
    float Vsettling;
} Data;

typedef struct
{
    float lowerBound;
    float upperBound;
    int newValue;
} Range;

float SettlingVelocity(int temp);
int RedefineTemp(float temp, Range ranges[]);
float Height(float alpha, float y);
float ShearVelocity(float Wspeed, float Precp);

int main()
{
    int i, j, k, m;
    float na, height, shearV;
    int one, two, three, four;
    int count = 0;
    float karman = 0.4;
    char line[100];
    short checkSum[200] = {0};
    float alphaOne, alphaTwo, beta, gamma, delta, alphaBetaOne, alphaBetaTwo;
    float flowDepthOne, flowDepthTwo, shearVelocity;
    float dnaOne, dnaTwo;
    float checkOne, checkTwo;

    Range ranges[] = {
        {7.5, 12.5, 10},
        {12.5, 17.5, 15},
        {17.5, 22.5, 20},
        {22.5, 27.5, 25},
        {27.5, 32.5, 30},
        {32.5, 37.5, 35},
    };

    Data data[200];


    //C0V370
    FILE *file = fopen("dataAllWeek_C0V310.csv", "r");

    if (file == NULL)
    {
        printf("Error opening file.\n");
        return 1;
    }

    fgets(line, 100, file);

    while(fgets(line, 100, file)!=NULL)
    {
      char *token;
      char *rest = line;

      // Extract and parse the index
      token = strtok_r(rest, ",", &rest);
      // We ignore the index (first column)

      // Extract and parse the date
      token = strtok_r(rest, ",", &rest);
      if (token != NULL) strcpy(data[count].Date, token);

      // Extract and parse DNALand
      token = strtok_r(rest, ",", &rest);
      if (token != NULL) data[count].DNALand = atof(token);

      // Extract and parse DNA06m
      token = strtok_r(rest, ",", &rest);
      if (token != NULL) data[count].DNA06m = atof(token);

      // Extract and parse DNA15m
      token = strtok_r(rest, ",", &rest);
      if (token != NULL) data[count].DNA15m = atof(token);

      // Extract and parse Temp
      token = strtok_r(rest, ",", &rest);
      if (token != NULL) data[count].Temp = atof(token);

      count++;
    }

    fclose(file);

    //redefine temperature
    for(i=0; i<count; i++)
    {
        data[i].Temp = RedefineTemp(data[i].Temp, ranges);
    }

    //calculate settling velocity
    for(i=0; i<count; i++)
    {
        data[i].Vsettling = SettlingVelocity(data[i].Temp);
    }

    printf("Date\t\tDNAland\t\tDNA0.6m\t\tDNA1.5m\t\tTemp\t\tVsettling\n");
    for(i=0; i<count; i++)
    {
        printf("%3d, %s\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.4f\n", i, data[i].Date, data[i].DNALand, data[i].DNA06m, data[i].DNA15m, data[i].Temp, data[i].Vsettling);
    }

    FILE *output_file = fopen("output.csv", "w");

    if (output_file == NULL) {
        printf("Error opening output file.\n");
        return 1;
    }

    fprintf(output_file, "date,na,flow depth,settling velocity,shear velocity,DNA0.6_cal,DNA0.6_theo,DNA1.5_cal,DNA1.5_theo\n");


    one = 0;
    na = 0.0012;
    beta = na / (1-na);
    height = 2.4;
    shearV = 0.081;

    for(i=0; i<count; i++) //date
    {
        two = 0;

        for(j=1; j<20; j++)
        //for(j=23; j<33; j++) //na
        {
            na = j / 10000.0;
            beta = na / (1-na);

            three = 0;

            for(k=15; k<=25; k++) //flow depth (h)
            {
                height = k / 10.0;

                flowDepthOne = 0.6/height;
                flowDepthTwo = 1.5/height;

                alphaOne = (1-flowDepthOne) / flowDepthOne;
                alphaTwo = (1-flowDepthTwo) / flowDepthTwo;

                alphaBetaOne = alphaOne * beta;
                alphaBetaTwo = alphaTwo * beta;

                four = 0;

                for(m=74; m<=84; m++) //shear velocity
                {
                    shearV = m/1000.0;

                    gamma = data[i].Vsettling / karman / shearV;

                    //gamma = data[i].Vsettling / karman / data[i].Vshear;

                    dnaOne = data[i].DNALand * pow(alphaBetaOne, gamma);
                    dnaTwo = data[i].DNALand * pow(alphaBetaTwo, gamma);

                    //checkOne = fabs(data[i].DNA06m - dnaOne);
                    //checkTwo = fabs(data[i].DNA15m - dnaTwo);

                    checkOne = fabs(data[i].DNA06m - dnaOne) / data[i].DNA06m;
                    checkTwo = fabs(data[i].DNA15m - dnaTwo) / data[i].DNA15m;

                    if(data[i].DNA06m != 0)
                    {
                        if(checkOne<0.1)
                        {
                            checkSum[one]++;

                            //printf("%f, %f\n", checkOne, checkTwo);
                            //printf("week: %d, na = %f, h = %f, u = %f \n\t 0.6c = %f, 0.6t = %f \t 1.5c = %f, 1.5t = %f \n\n", i, na, height, shearV, dnaOne, data[i].DNA06m, dnaTwo, data[i].DNA15m);
                            fprintf(output_file, "%s,%f,%f,%f,%f,%f,%f,%f,%f\n", data[i].Date, na, height, data[i].Vsettling, shearV, dnaOne, data[i].DNA06m, dnaTwo, data[i].DNA15m);
                        }
                    }

                    //if(m%10 == 0)
                    {
                        four++;
                    }
                }

                //if(k%5 == 0)
                {
                    three++;
                }
            }
            two++;
        }
        one++;
    }


/*
    double deltaOne, deltaTwo;
    double heightOne, heightTwo;
    double shearVOne, shearVTwo;

    flowDepthOne = 0.6/height;
    flowDepthTwo = 1.5/height;

    alphaOne = (1-flowDepthOne) / flowDepthOne;
    alphaTwo = (1-flowDepthTwo) / flowDepthTwo;

    alphaBetaOne = alphaOne * beta;
    alphaBetaTwo = alphaTwo * beta;

    alphaBetaOne = log(alphaBetaOne);
    alphaBetaTwo = log(alphaBetaTwo);

    for(i=21; i<40; i++) //week
    {
        deltaOne = data[i].DNA06m / data[i].DNALand;
        deltaTwo = data[i].DNA15m / data[i].DNALand;

        deltaOne = log(deltaOne);
        deltaTwo = log(deltaTwo);

        gamma = data[i].Vsettling / karman;

        shearVOne = gamma * alphaBetaOne / deltaOne;
        shearVTwo = gamma * alphaBetaTwo / deltaTwo;

        //gamma =  karman * data[i].Vshear / data[i].Vsettling;

        //deltaOne = pow(deltaOne, gamma);
        //deltaTwo = pow(deltaTwo, gamma);

        //alphaOne = deltaOne / beta;
        //alphaTwo = deltaTwo / beta;

        //heightOne = Height(alphaOne, 0.6);
        //heightTwo = Height(alphaTwo, 1.5);

        printf("%d: %f\t%f\n", i, shearVOne, shearVTwo);
    }
*/

    fclose(output_file);

    printf("CSV file 'output.csv' has been created successfully.\n");

    printf("\n\n%d\t%d\t%d\t%d\n", one, two, three, four);

    one = 0;

    for(i=0; i<count; i++)
    {
        if(checkSum[i] > 0)
        {
            printf("%d, %d\n", checkSum[i], i);
            one++;
        }
    }
    printf("\n%d\n", one);

    return 0;
}

int RedefineTemp(float temp, Range ranges[])
{
    int i;

    for(i=0; i<5; i++)
    {
        if(temp>=ranges[i].lowerBound && temp<ranges[i].upperBound)
        {
            return ranges[i].newValue;
        }
    }
    return temp;
}

float SettlingVelocity(int temp)
{
    if(temp == 10)
    {
        return 0.0051151;
    }else if(temp == 15)
    {
        return 0.0050471;
    }else if(temp == 20)
    {
        return 0.0049835;
    }else if(temp == 25)
    {
        return 0.0049189;
    }else if(temp == 30)
    {
        return 0.0048585;
    }else if(temp == 35)
    {
        return 0.0047996;
    }

    return 0;
}

float Height(float alpha, float y)
{
    float n;

    n = 1 / (alpha+1);

    return y / n;
}
