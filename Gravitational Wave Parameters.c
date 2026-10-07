#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
 
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define N 2000 //Number of data points
#define M 200000 //Number of random walk points
#define B 200 //Number of bins (for plotting)

// I have used AI to generate the normal distribution, since proper random number generators do not seem to exist by default in C
static uint64_t s[4];

static uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline uint64_t rotl(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t next_u64(void) { //AI-generated xoshiro256** algorithm
    uint64_t result = rotl(s[1] * 5, 7) * 9;
    uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1];
    s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return result;
}

void rng_seed(uint64_t seed) {
    for (int i = 0; i < 4; i++) s[i] = splitmix64(&seed);
}

double rng_uniform(void) {/* uniform in (0,1), 53-bit */
    return ((next_u64() >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

double rng_normal(double mu, double sigma) {/* Box-Muller, caches spare */
    static int have = 0;
    static double spare;
    if (have) { have = 0; return mu + sigma * spare; }
    double u1 = rng_uniform(), u2 = rng_uniform();
    double r = sqrt(-2.0 * log(u1)), th = 2.0 * M_PI * u2;
    spare = r * sin(th);
    have = 1;
    return mu + sigma * r * cos(th);
}
//This is till where the AI generated part is

double h(double a, double b, double g, double t){ //Gravitational wave function h(t) with the parameters
    return a*exp(t)*(1-tanh(2*(t-b)))*sin(g*t);
}

void make_pdf(double *var, double lo, double hi, double *varc, double *pdfvar){ //Sorting the values of alpha, beta and gamma we got from the random walk into the respective bins
    long cnt[B]={0}, tot=0;
    double w=(hi-lo)/B; //Size of each bin
    for(int j=70000; j<M; j++){ //70000 is the burn-in
        int bin=(int)((var[j]-lo)/w);
        if(bin>=0 && bin<B){cnt[bin]++; tot++;}
    }
    for(int i=0; i<B; i++){
        varc[i]=lo+(i+0.5)*w; //Centre points of each bin
        pdfvar[i]=cnt[i]/(tot*w); //Normalising the values
    }
}

void write_pdf(const char *fname, double *varc, double *pdfvar, int n){ //Setup to plot the graph
    FILE *out=fopen(fname, "w"); if(!out){perror("fopen"); return;}
    for (int i=0; i<n; i++) fprintf(out, "%f %f\n", varc[i], pdfvar[i]);
    fclose(out);
}

void plot_pdf(const char *fname, const char *name){ //Function to actually plot the graph
    FILE *gp=popen("gnuplot -persist", "w"); if(!gp){perror("popen"); return;}
    fprintf(gp, "set title 'PDF of %s'\n", name);
    fprintf(gp, "set xlabel '%s'\n", name);
    fprintf(gp, "set ylabel 'P(%s|data)'\n", name);
    fprintf(gp, "set style fill solid 0.5\n");
    fprintf(gp, "plot '%s' with boxes notitle\n", fname);
    pclose(gp);
}

static double tval[N], hval[N]; //data points (given)
static double a[M], b[M], g[M]; //arrrays to store alpha, beta and gamma

double MHAlgo(double al, double be, double ga){
    double Y=0.0;
    for(int k=0; k<N; k++){
        double temp=(hval[k]-h(al,be,ga,tval[k]))/0.2;
        Y+=temp*temp;
    }
    return -Y;
}

int main(){
    rng_seed((uint64_t)time(NULL) ^ ((uint64_t)clock() << 32)); //initialising the RNG (also AI-generated)

    //Storing the values of time and h from the given data:
    FILE *f = fopen("gw_data.csv", "r"); if (!f){perror("fopen"); return 1;}
    char line[256]; fgets(line, sizeof line, f);
    size_t n=0;
    while(n<N && fgets(line, sizeof line, f)){
        if (sscanf(line, "%lf,,%lf", &tval[n], &hval[n])==2) n++;
    }
    fclose(f);

    //Executing the M-H algorithm:
    double ca=1.9, cb=4.6, cg=11.4636242; //Randomly selected starting points for alpha, beta and gamma respectively
    double Y=MHAlgo(ca, cb, cg);
    for(int j=0; j<M; j++){
        double pa = ca + rng_normal(0.0, 0.000079); //The value of mu must be 0, and the values of sigma were taken after trial and error
        double pb = cb + rng_normal(0.0, 0.00039);
        double pg = cg + rng_normal(0.0, 0.00019);
        if(pa>0 && pa<2 && pb>1 && pb<10 && pg>1 && pg<20){
            double Yn=MHAlgo(pa, pb, pg);
            if(rng_uniform()<exp(Yn-Y)){ca=pa; cb=pb; cg=pg; Y=Yn;} //Takes new values if the step is accepted
        }
        a[j]=ca; b[j]=cb; g[j]=cg; //Stores the values at each step of the random walk
    }
    
    //Plotting the PDFs:
    double ac[B], bc[B], gc[B], pdfa[B], pdfb[B], pdfg[B]; //Initialising the centre points and PDFs of each parameter
    make_pdf(a, 0.0, 2.0, ac, pdfa);
    make_pdf(b, 1.0, 10.0, bc, pdfb);
    make_pdf(g, 1.0, 20.0, gc, pdfg);
    write_pdf("alpha_pdf.dat", ac, pdfa, B);
    write_pdf("beta_pdf.dat", bc, pdfb, B);
    write_pdf("gamma_pdf.dat", gc, pdfg, B);
    plot_pdf("alpha_pdf.dat", "alpha");
    plot_pdf("beta_pdf.dat", "beta");
    plot_pdf("gamma_pdf.dat", "gamma");
   
    return 0;
}