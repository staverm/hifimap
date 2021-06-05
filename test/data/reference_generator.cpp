#include <bits/stdc++.h>
#include <sys/time.h>

void generate_reference(int REFERENCE_SIZE, int SEED, const std::map<int,int>& add_seeds) {
    std::default_random_engine generator;
    generator.seed(SEED);
    std::default_random_engine generator_err;
    generator_err.seed(SEED*SEED);

    std::uniform_real_distribution<float> fl;

    std::uniform_int_distribution<int> bases_distribution(0, 3);
    std::vector<char> bases = {'A', 'C', 'G', 'T'};
    printf(">ARTIFICIAL_REF\n");
    int print_new_line = 0;
    for (int i = 0; i < REFERENCE_SIZE; i++) {
        print_new_line++;
        if (add_seeds.count(i)){
            generator.seed(SEED+add_seeds.at(i));
        }
        int ind = bases_distribution(generator);
        if (fl(generator_err)<0.01) ind = (ind+bases_distribution(generator_err))%4;
        std::cout << bases[ind];
        if(print_new_line==80) {
            std::cout << std::endl;
            print_new_line=0;
        }
    }   
    printf("\n");
}

int main(int argc, char** argv) {
    if (argc < 2 + 1) {
        std::cout << "wrong args, usage: REFERENCE_SIZE SEED" << std::endl;
        exit(1);
    }
    if (argc == 3){
        generate_reference(atoi(argv[1]), atoi(argv[2]), std::map<int,int>{});
    }
    else{
        int map_size = atoi(argv[3]);
        std::map<int,int> m;
        for (int i=0;i<map_size;i++){
            m[atoi(argv[4+2*i])]=atoi(argv[4+2*i+1]);
        }
        generate_reference(atoi(argv[1]), atoi(argv[2]), m);
    }
}