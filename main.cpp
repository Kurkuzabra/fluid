

#include <bits/stdc++.h>
#include "liq.cpp"

using namespace std;


#define FLOAT float
#define DOUBLE double
#define FIXED(n, k) Fixed<n, k>
#define FAST_FIXED(n, k) Fixed<n, k, true>
#define S(n, k) pair<uint32_t, uint32_t>(n, k)

#ifndef TYPES
#error "Not defined any types"
#endif

#ifndef SIZES
#error "Not defined any sizes"
#endif



// constexpr size_t N = 36, M = 84;
// constexpr size_t N = 14, M = 5;
constexpr size_t T = 1000;

// char field[N][M + 1] = {
//     "#####",
//     "#.  #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#...#",
//     "#####",
//     "#   #",
//     "#   #",
//     "#   #",
//     "#####",
// };

// char field[N][M + 1] = {
//     "####################################################################################",
//     "#                                                                                  #",
//     "#                                                                                  #",
//     "#                                                                                  #",
//     "#                                                                                  #",
//     "#                                                                                  #",
//     "#                                       .........                                  #",
//     "#..............#            #           .........                                  #",
//     "#..............#            #           .........                                  #",
//     "#..............#            #           .........                                  #",
//     "#..............#            #                                                      #",
//     "#..............#            #                                                      #",
//     "#..............#            #                                                      #",
//     "#..............#            #                                                      #",
//     "#..............#............#                                                      #",
//     "#..............#............#                                                      #",
//     "#..............#............#                                                      #",
//     "#..............#............#                                                      #",
//     "#..............#............#                                                      #",
//     "#..............#............#                                                      #",
//     "#..............#............#                                                      #",
//     "#..............#............#                                                      #",
//     "#..............#............################                     #                 #",
//     "#...........................#....................................#                 #",
//     "#...........................#....................................#                 #",
//     "#...........................#....................................#                 #",
//     "##################################################################                 #",
//     "#                                                                                  #",
//     "#                                                                                  #",
//     "#                                                                                  #",
//     "#                                                                                  #",
//     "#                                                                                  #",
//     "#                                                                                  #",
//     "#                                                                                  #",
//     "#                                                                                  #",
//     "####################################################################################",
// };

template<typename arg1, typename arg2, typename arg3, pair<uint32_t, uint32_t> sz>
constexpr void declare_class()
{
    liq::FieldEmulator<arg1, arg2, arg3, sz.first, sz.second>::declare();
}

template<typename arg1, typename arg2, typename arg3, pair<uint32_t, uint32_t>... Args>
constexpr void declare_fields()
{
    ((void)declare_class<arg1, arg2, arg3, Args>(), ...);
}

template<typename arg1, typename arg2, typename... Args>
constexpr void decl_types_1()
{
    ((void)declare_fields<arg1, arg2, Args, SIZES>(), ...);
}

template<typename arg1, typename... Args>
constexpr void decl_types_2()
{
    ((void)decl_types_1<arg1, Args, Args...>(), ...);
}

template<typename... Args>
constexpr bool decl_types_3()
{
    ((void)decl_types_2<Args, Args...>(), ...);
    return true;
}

constexpr bool compiled = decl_types_3<TYPES>();

int main()
{
    string line;
    vector<string> field;
    std::ifstream myfile("field.txt");

    while (std::getline(myfile, line))
    {
        field.push_back(line);
    }

    liq::DynamicFieldEmulator<Fixed<32, 16>, Fixed<64, 16>, Fixed<16, 16>> emul = 
        liq::DynamicFieldEmulator<Fixed<32, 16>, Fixed<64, 16>, Fixed<16, 16>>(14, 10, field);
      
    liq::FieldEmulator<Fixed<32, 16>, Fixed<64, 16>, Fixed<16, 16>, 14, 10> emul_ = 
        liq::FieldEmulator<Fixed<32, 16>, Fixed<64, 16>, Fixed<16, 16>, 14, 10>(field);

    emul.set_vel(' ', 0.01);
    emul.set_vel('.', 1000);
    emul.exec(T, 0.1);    

    emul_.set_vel(' ', 0.01);
    emul_.set_vel('.', 1000);
    emul_.exec(T, 0.1);  

}