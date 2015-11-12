#include <stdlib.h>     /* strtoul */
#include "include/para.h"
#include "getopt.h"

// Totally 12 global varibles (parameters) in the following
char   	INPUT_PATH[500];
char   	INPUT_FILE[500];
char   	OUTPUT_PATH[500];
char   	OUTPUT_FILE[500];

double KRYLOV_TOL;
unsigned long KRYLOV_M;
unsigned long KRYLOV_ITRACE;

unsigned long KLIM;
unsigned long MAX_THREADS_PER_BLOCK;
unsigned long MAX_GRID_SIZE_1;
unsigned long MAX_GRID_SIZE_2;
unsigned long MAX_GRID_SIZE_3;

void print_help();

static struct option long_options[] = 
{
        {"iPATH",     required_argument, 0,  0 },//#01
        {"iFILE",     required_argument, 0,  0 },//#02
        {"oPATH",     required_argument, 0,  0 },//#03
        {"oFILE",     required_argument, 0,  0 },//#04
        {"TOL",       required_argument, 0,  0 },//#05
        {"KRYLOV_M",  required_argument, 0,  0 },//#06
        {"ITRACE",    required_argument, 0,  0 },//#07
        {"KLIM",      required_argument, 0,  0 },//#08
        {"THREADS",   required_argument, 0,  0 },//#09
        {"GRID1",     required_argument, 0,  0 },//#10
        {"GRID2",     required_argument, 0,  0 },//#11
        {"GRID3",     required_argument, 0,  0 },//#12
        {"help",  no_argument,       0,  0 }
};

void set_parameter_value(char * name, char * val)
{
    //I/O /////////////////////////////
        if(strcmp("iPATH",name) == 0)
        {
                strcpy(INPUT_PATH, val);
        }
        if(strcmp("iFILE",name) == 0)
        {
                strcpy(INPUT_FILE, val);
        }
        if(strcmp("oPATH",name) == 0)
        {
                strcpy(OUTPUT_PATH, val);
        }
        if(strcmp("oFILE",name) == 0)
        {
                strcpy(OUTPUT_FILE, val);
        }
        
    //Krylov space /////////////////////////////
        if(strcmp("TOL",name) == 0)//#01
        {
                KRYLOV_TOL = atof(val);
        }
        if(strcmp("KRYLOV_M",name) == 0)//#01
        {
                KRYLOV_M = strtoul(val, NULL, 0);
        }
        if(strcmp("ITRACE",name) == 0)//#01
        {
                KRYLOV_ITRACE = strtoul(val, NULL, 0);
        }
        
    //CUDAe /////////////////////////////
        if(strcmp("KLIM",name) == 0)//#01
        {
                KLIM = strtoul(val, NULL, 0);
        }
        if(strcmp("THREADS",name) == 0)//#01
        {
                MAX_THREADS_PER_BLOCK= strtoul(val, NULL, 0);
        }
        if(strcmp("GRID1",name) == 0)//#01
        {
                MAX_GRID_SIZE_1= strtoul(val, NULL, 0);
        }
        if(strcmp("GRID2",name) == 0)//#01
        {
                MAX_GRID_SIZE_2= strtoul(val, NULL, 0);
        }
        if(strcmp("GRID3",name) == 0)//#01
        {
                MAX_GRID_SIZE_3= strtoul(val, NULL, 0);
        }
        if(strcmp("help",name) == 0)//#01
        {
            print_help();
        }
}

void Load_Config()
{
        char *cfg_path;
        char cfg[200], name[100], val[200];

        cfg_path= new char [100];
        cfg_path= getenv("EVO_CFG_PATH");
        strcpy(cfg, cfg_path);
        strcat(cfg, "/Config.cfg");

        ifstream config(cfg);
        if (!config)
        {
                cout << "***** ERROR ***** " << endl;
                cout << "Cannot load config. file: " << cfg << endl;
                assert(0);
        }

        int nline, i;
        string line;
        char * cstr, *token, *cp;
        const char delimiters[] = "= \t";

        config >> nline;
        for (i = 0; i < nline+1; i++)
        {
                getline(config, line);
                if (i > 0)
                {
                        cstr = new char [line.size() + 1];
                        strcpy(cstr, line.c_str());

                        cp = strdup(cstr);

                        token = strtok(cp, delimiters);
                        strcpy(name, token);

                        token = strtok(NULL, delimiters);
                        strcpy(val, token);

                        delete [] cstr;

                        set_parameter_value(name, val);
                }
        }
        config.close();
}

void ParameterResolve(int argc, char ** argv)
{
        int c;
        int digit_optind = 0;
        char name[100];

        char opt_str[500];
        strcpy(opt_str, "");
        //////////////////////////////////////////////////////
        //load config
        Load_Config();

        optind=1;
        while (1) {
                int this_option_optind = optind ? optind : 1;
                int option_index = 0;

                c = getopt_long(argc, argv, "h", long_options, &option_index);
                if (c == -1) break;

                switch (c) {
                        case 0:
//                                printf("option %s", long_options[option_index].name);
//                                if (optarg)
//                                        printf(" with arg %s", optarg);
//                                printf("\n");
                                strcpy(name, long_options[option_index].name);
                                
                                set_parameter_value(name, optarg);

                                if (strcmp("oFILE", name)!=0)
                                {
                                    strcat(opt_str, "_");
                                    strcat(opt_str, name);
                                    strcat(opt_str, "_");
                                    strcat(opt_str, optarg);
                                }
                                break;

                        case 'h':
                                print_help();
                        case '?':
                                break;

                        default:
                                printf("?? getopt returned character code 0%o ??\n", c);
                }

        }
        strcat(OUTPUT_FILE, opt_str);
        strcat(OUTPUT_FILE, ".dat");
}

void PrintParameters()
{
        cout << "========================================================="<< endl;
        cout << "Input Path=" << "\t" << INPUT_PATH << endl;
        cout <<	"Input File=" << "\t" << INPUT_FILE << endl;
        cout <<	"Output Path=" << "\t" << OUTPUT_PATH << endl;
        cout <<	"Output File=" << "\t" << OUTPUT_FILE << endl;
        
        cout << "tol=" << "\t" << KRYLOV_TOL << endl;
        cout << "m =" << "\t" << KRYLOV_M << endl;
        cout << "itrace=" << "\t" << KRYLOV_ITRACE << endl;
        
        cout << "klim=" << "\t" << KLIM << endl;
        cout << "threads=" << "\t" << MAX_THREADS_PER_BLOCK << endl;
        cout << "grid1=" << "\t" << MAX_GRID_SIZE_1 << endl;
        cout << "grid2=" << "\t" << MAX_GRID_SIZE_2 << endl;
        cout << "grid3=" << "\t" << MAX_GRID_SIZE_3 << endl;
        cout << "========================================================="<< endl;
        cout << endl;
}


void print_help()
{
    cout << "Please use the following options:" << endl;
    cout << endl;
    cout << "--iPATH=path" << "\t\t" << "specifies the INPUT_PATH" << endl;
    cout << "--iFILE=file" << "\t\t" << "specifies the INPUT_FILE" << endl;
    cout << "--oPATH=path" << "\t\t" << "specifies the OUTPUT_PATH" << endl;
    cout << "--oFILE=file" << "\t\t" << "specifies the OUTPUT_FILE" << endl;

    cout << "--TOL=val (doulbe)" << "\t" << "specifies the tolerance value of the Krylov space method."  << endl;
    cout << "--KRYLOV_M=val (int)" << "\t" << "specifies the dimension of the Krylov space." << endl;
    cout << "--ITRACE=val (bool)" << "\t" << "flag for printing information of Krylov space method." << endl;

    cout << "--KLIM=val (int)" << "\t" << "specifies the K-limite of GPU implementation." << endl;
    cout << "--THREADS=val (int)" << "\t" << "specifies the MAX_THREADS_PER_BLOCK." << endl;
    cout << "--GRID1=val (int)" << "\t" << "specifies MAX_GRID_SIZE_1." << endl;
    cout << "--GRID2=val (int)" << "\t" << "specifies MAX_GRID_SIZE_2." << endl;
    cout << "--GRID3=val (int)" << "\t" << "specifies MAX_GRID_SIZE_3." << endl;
    cout << endl;
    exit(EXIT_SUCCESS);
}
