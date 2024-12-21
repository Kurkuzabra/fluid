#include "fixed.cpp"
#include <bits/stdc++.h>

using namespace std;

constexpr size_t N_ = 14, M_ = 10;
// char field[N_][M_ + 1] = {
//     "##########",
//     "#        #",
//     "# .....  #",
//     "# .....  #",
//     "# .....  #",
//     "#        #",
//     "#        #",
//     "#        #",
//     "#        #",
//     "#        #",
//     "#        #",
//     "#        #",
//     "#        #",
//     "##########",
// };

namespace liq
{
    constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};
    static constexpr Fixed<32, 16> inf = Fixed<32, 16>::from_raw(std::numeric_limits<int32_t>::max());
    static constexpr Fixed<32, 16> eps = Fixed<32, 16>::from_raw(deltas.size());
    mt19937 rnd(1337);

    template<typename T, uint32_t N, uint32_t M>
    struct VectorField
    {
        array<T, deltas.size()> v[N][M];
        T &add(int x, int y, int dx, int dy, T dv) {
            return get(x, y, dx, dy) += dv;
        }

        T &get(int x, int y, int dx, int dy) {
            size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
            assert(i < deltas.size());
            return v[x][y][i];
        }
    };

    template<typename T>
    struct DynamicVectorField
    {
        vector<vector<array<T, deltas.size()>>> v;
        DynamicVectorField(uint32_t N, uint32_t M)
        {
            v.resize(N);
            for (uint32_t i = 0; i < N; i++)
            {
                v[i].resize(M);
            }
        }
        T &add(int x, int y, int dx, int dy, T dv) {
            return get(x, y, dx, dy) += dv;
        }

        T &get(int x, int y, int dx, int dy) {
            size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
            assert(i < deltas.size());
            return v[x][y][i];
        }
    };

    template<typename PType, typename VType>
    struct ParticleParams
    {
        char type;
        PType cur_p;
        array<VType, deltas.size()> v;
    };

    template<typename PType, typename VType, typename VFType, uint32_t N, uint32_t M>
    class FieldEmulator
    {
        PType rho[256];

        PType p[N][M]{}, old_p[N][M];
        int dirs[N][M]{};

        VectorField<VType, N, M> velocity{};
        VectorField<VFType, N, M> velocity_flow{};
        vector<string> field;
        int last_use[N][M]{};
        int UT = 0; 

        void swap_with(ParticleParams<PType, VType>& params, int x, int y) {
            swap(field[x][y], params.type);
            swap(p[x][y], params.cur_p);
            swap(velocity.v[x][y], params.v);
        }


    public:

        static constexpr void declare() {}

        // FieldEmulator(char field_[N][M + 1]) : field(field_)
        FieldEmulator(vector<string> field_) : field(field_)
        {
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    // cout << field[x][y] << " ";
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy] : liq::deltas) {
                        dirs[x][y] += (field[x + dx][y + dy] != '#');
                    }
                    // cout << dirs[x][y] << "   ";
                }
                // cout << endl;
            }
            // cout << "created\n";
        }

        void set_vel(char let, double k)
        {
            rho[let] = k;
        }

            tuple<VType, bool, pair<int, int>> propagate_flow(int x, int y, VType lim) {
                // cout << lim;
                last_use[x][y] = UT - 1;
                VType ret = 0;
                for (auto [dx, dy] : deltas) {
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
                        VType cap = velocity.get(x, y, dx, dy);
                        VType flow = velocity_flow.get(x, y, dx, dy);
                        if (flow == cap) {
                            continue;
                        }
                        // assert(v >= velocity_flow.get(x, y, dx, dy));
                        auto vp = min(lim, cap - flow);
                        if (last_use[nx][ny] == UT - 1) {
                            velocity_flow.add(x, y, dx, dy, vp);
                            last_use[x][y] = UT;
                            // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
                            // cout << nx << " " << ny << endl;
                            return {vp, 1, {nx, ny}};
                        }
                        auto [t, prop, end] = propagate_flow(nx, ny, vp);
                        ret += t;
                        if (prop) {
                            velocity_flow.add(x, y, dx, dy, t);
                            last_use[x][y] = UT;
                            // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t << " / " << lim << "\n";
                            // cout << end.first << " " << end.second << endl;
                            return {t, prop && end != pair(x, y), end};
                        }
                    }
                }
                last_use[x][y] = UT;
                return {ret, 0, {0, 0}};
            }

            void propagate_stop(int x, int y, bool force = false) {
                if (!force) {
                    bool stop = true;
                    for (auto [dx, dy] : deltas) {
                        int nx = x + dx, ny = y + dy;
                        if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > 0) {
                            stop = false;
                            break;
                        }
                    }
                    if (!stop) {
                        return;
                    }
                }
                last_use[x][y] = UT;
                for (auto [dx, dy] : deltas) {
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > 0) {
                        continue;
                    }
                    propagate_stop(nx, ny);
                }
            }

            VType move_prob(int x, int y) {
                VType sum = 0;
                for (size_t i = 0; i < deltas.size(); ++i) {
                    auto [dx, dy] = deltas[i];
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                        continue;
                    }
                    auto v = velocity.get(x, y, dx, dy);
                    if (v < 0) {
                        continue;
                    }
                    sum += v;
                }
                return sum;
            }

            bool propagate_move(int x, int y, bool is_first) {
                last_use[x][y] = UT - is_first;
                bool ret = false;
                int nx = -1, ny = -1;
                do {
                    std::array<VType, deltas.size()> tres;
                    VType sum = 0;
                    for (size_t i = 0; i < deltas.size(); ++i) {
                        auto [dx, dy] = deltas[i];
                        int nx = x + dx, ny = y + dy;
                        if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                            tres[i] = sum;
                            continue;
                        }
                        auto v = velocity.get(x, y, dx, dy);
                        if (v < 0) {
                            tres[i] = sum;
                            continue;
                        }
                        sum += v;
                        tres[i] = sum;
                    }

                    if (sum == 0) {
                        break;
                    }

                    VType p = VType::random01(rnd) * sum;
                    size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

                    auto [dx, dy] = deltas[d];
                    nx = x + dx;
                    ny = y + dy;
                    assert(velocity.get(x, y, dx, dy) > 0 && field[nx][ny] != '#' && last_use[nx][ny] < UT);

                    ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
                } while (!ret);
                last_use[x][y] = UT;
                for (size_t i = 0; i < deltas.size(); ++i) {
                    auto [dx, dy] = deltas[i];
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < 0) {
                        propagate_stop(nx, ny);
                    }
                }
                if (ret) {
                    if (!is_first) {
                        ParticleParams<PType, VType> pp{};
                        swap_with(pp, x, y);
                        swap_with(pp, nx, ny);
                        swap_with(pp, x, y);
                    }
                }
                return ret;
            }
        void exec(uint32_t Time, VType g = 0.1)
        {
            cout << "exec\n";
            for (size_t i = 0; i < Time; ++i) {
                
                PType total_delta_p = 0;
                // Apply external forces
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] == '#')
                            continue;
                        if (field[x + 1][y] != '#')
                            velocity.add(x, y, 1, 0, g);
                    }
                }

                // Apply forces from p
                memcpy(old_p, p, sizeof(p));
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] == '#')
                            continue;
                        for (auto [dx, dy] : deltas) {
                            int nx = x + dx, ny = y + dy;
                            if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
                                PType delta_p = old_p[x][y] - old_p[nx][ny];
                                PType force = delta_p;
                                VType &contr = velocity.get(nx, ny, -dx, -dy);
                                if (contr * rho[(int) field[nx][ny]] >= force) {
                                    contr -= force / rho[(int) field[nx][ny]];
                                    continue;
                                }
                                force -= contr * rho[(int) field[nx][ny]];
                                contr = 0;
                                velocity.add(x, y, dx, dy, force / rho[(int) field[x][y]]);
                                p[x][y] -= force / dirs[x][y];
                                total_delta_p -= force / dirs[x][y];
                            }
                        }
                        // cout << "deltas\n";
                    }
                }

                // Make flow from velocities
                velocity_flow = {};
                bool prop = false;
                do {
                    UT += 2;
                    prop = 0;
                    for (size_t x = 0; x < N; ++x) {
                        for (size_t y = 0; y < M; ++y) {
                            if (field[x][y] != '#' && last_use[x][y] != UT) {
                                auto [t, local_prop, _] = propagate_flow(x, y, 1);
                                if (t > 0) {
                                    prop = 1;
                                }
                            }
                        }
                    }
                } while (prop);

                // Recalculate p with kinetic energy
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] == '#')
                            continue;
                        for (auto [dx, dy] : deltas) {
                            VType old_v = velocity.get(x, y, dx, dy);
                            VType new_v = velocity_flow.get(x, y, dx, dy);
                            if (old_v > 0) {
                                // cout << new_v.v << " " << old_v.v << endl;
                                assert(new_v <= old_v);
                                velocity.get(x, y, dx, dy) = new_v;
                                PType force = (old_v - new_v) * rho[(int) field[x][y]];
                                if (field[x][y] == '.')
                                    force *= 0.8;
                                if (field[x + dx][y + dy] == '#') {
                                    // cout << "# ";
                                    p[x][y] += force / dirs[x][y];
                                    total_delta_p += force / dirs[x][y];
                                } else {
                                    // cout << "!# ";
                                    p[x + dx][y + dy] += force / dirs[x + dx][y + dy];
                                    total_delta_p += force / dirs[x + dx][y + dy];
                                }
                            }
                            // cout << "recalc\n";
                        }
                    }
                }

                UT += 2;
                prop = false;
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] != '#' && last_use[x][y] != UT) {
                            if (VType::random01(rnd) < move_prob(x, y)) {
                                prop = true;
                                propagate_move(x, y, true);
                            } else {
                                propagate_stop(x, y, true);
                            }
                        }
                    }
                }

                // cout << "Tick " << i << ":\n";
                // for (size_t x = 0; x < N; ++x) {
                //     for (size_t y = 0; y < M; ++y) {
                //         cout << p[x][y] << " ";
                //     }
                //     cout << endl;
                // }
                // cout << endl;
                if (prop) {
                    cout << "Tick " << i << ":\n";
                    for (size_t x = 0; x < N; ++x) {
                        cout << field[x] << "\n";
                    }
                }
            }
        }

    };

    template<typename PType, typename VType, typename VFType>
    class DynamicFieldEmulator
    {
        uint32_t N, M;
        PType rho[256];

        vector<vector<PType>> p, old_p;
        vector<vector<int>> dirs;

        DynamicVectorField<VType> velocity;
        DynamicVectorField<VFType> velocity_flow;
        vector<string> field;
        vector<vector<int>> last_use;
        int UT = 0; 

        void swap_with(ParticleParams<PType, VType>& params, int x, int y) {
            swap(field[x][y], params.type);
            swap(p[x][y], params.cur_p);
            swap(velocity.v[x][y], params.v);
        }


    public:

        static constexpr void declare() {}
        DynamicFieldEmulator(uint32_t N__, uint32_t M__, vector<string> field_) : field(field_), velocity(DynamicVectorField<VType>(N__, M__)), velocity_flow(DynamicVectorField<VFType>(N__, M__))
        {
            N = N__;
            M = M__;
            p.resize(N);
            old_p.resize(N);
            dirs.resize(N);
            last_use.resize(N);
            for (uint32_t i = 0; i < N; i++)
            {
                p[i].resize(M);
                old_p[i].resize(M);
                dirs[i].resize(M);
                last_use[i].resize(M);
            }
            // velocity = DynamicVectorField<VType>(N, M);
            // velocity_flow = DynamicVectorField<VFType>(N, M);


            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    // cout << field[x][y] << " ";
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy] : liq::deltas) {
                        dirs[x][y] += (field[x + dx][y + dy] != '#');
                    }
                    // cout << dirs[x][y] << "   ";
                }
                // cout << endl;
            }
            // cout << "created\n";
        }

        void set_vel(char let, double k)
        {
            rho[let] = k;
        }

            tuple<VType, bool, pair<int, int>> propagate_flow(int x, int y, VType lim) {
                // cout << lim;
                last_use[x][y] = UT - 1;
                VType ret = 0;
                for (auto [dx, dy] : deltas) {
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
                        VType cap = velocity.get(x, y, dx, dy);
                        VType flow = velocity_flow.get(x, y, dx, dy);
                        if (flow == cap) {
                            continue;
                        }
                        // assert(v >= velocity_flow.get(x, y, dx, dy));
                        auto vp = min(lim, cap - flow);
                        if (last_use[nx][ny] == UT - 1) {
                            velocity_flow.add(x, y, dx, dy, vp);
                            last_use[x][y] = UT;
                            // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
                            // cout << nx << " " << ny << endl;
                            return {vp, 1, {nx, ny}};
                        }
                        auto [t, prop, end] = propagate_flow(nx, ny, vp);
                        ret += t;
                        if (prop) {
                            velocity_flow.add(x, y, dx, dy, t);
                            last_use[x][y] = UT;
                            // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t << " / " << lim << "\n";
                            // cout << end.first << " " << end.second << endl;
                            return {t, prop && end != pair(x, y), end};
                        }
                    }
                }
                last_use[x][y] = UT;
                return {ret, 0, {0, 0}};
            }

            void propagate_stop(int x, int y, bool force = false) {
                if (!force) {
                    bool stop = true;
                    for (auto [dx, dy] : deltas) {
                        int nx = x + dx, ny = y + dy;
                        if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > 0) {
                            stop = false;
                            break;
                        }
                    }
                    if (!stop) {
                        return;
                    }
                }
                last_use[x][y] = UT;
                for (auto [dx, dy] : deltas) {
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > 0) {
                        continue;
                    }
                    propagate_stop(nx, ny);
                }
            }

            VType move_prob(int x, int y) {
                VType sum = 0;
                for (size_t i = 0; i < deltas.size(); ++i) {
                    auto [dx, dy] = deltas[i];
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                        continue;
                    }
                    auto v = velocity.get(x, y, dx, dy);
                    if (v < 0) {
                        continue;
                    }
                    sum += v;
                }
                return sum;
            }

            bool propagate_move(int x, int y, bool is_first) {
                last_use[x][y] = UT - is_first;
                bool ret = false;
                int nx = -1, ny = -1;
                do {
                    std::array<VType, deltas.size()> tres;
                    VType sum = 0;
                    for (size_t i = 0; i < deltas.size(); ++i) {
                        auto [dx, dy] = deltas[i];
                        int nx = x + dx, ny = y + dy;
                        if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                            tres[i] = sum;
                            continue;
                        }
                        auto v = velocity.get(x, y, dx, dy);
                        if (v < 0) {
                            tres[i] = sum;
                            continue;
                        }
                        sum += v;
                        tres[i] = sum;
                    }

                    if (sum == 0) {
                        break;
                    }

                    VType p = VType::random01(rnd) * sum;
                    size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

                    auto [dx, dy] = deltas[d];
                    nx = x + dx;
                    ny = y + dy;
                    assert(velocity.get(x, y, dx, dy) > 0 && field[nx][ny] != '#' && last_use[nx][ny] < UT);

                    ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
                } while (!ret);
                last_use[x][y] = UT;
                for (size_t i = 0; i < deltas.size(); ++i) {
                    auto [dx, dy] = deltas[i];
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < 0) {
                        propagate_stop(nx, ny);
                    }
                }
                if (ret) {
                    if (!is_first) {
                        ParticleParams<PType, VType> pp{};
                        swap_with(pp, x, y);
                        swap_with(pp, nx, ny);
                        swap_with(pp, x, y);
                    }
                }
                return ret;
            }
        void exec(uint32_t Time, VType g = 0.1)
        {
            cout << "exec\n";
            for (size_t i = 0; i < Time; ++i) {
                
                PType total_delta_p = 0;
                // Apply external forces
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] == '#')
                            continue;
                        if (field[x + 1][y] != '#')
                            velocity.add(x, y, 1, 0, g);
                    }
                }

                // Apply forces from p
                old_p = p;
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] == '#')
                            continue;
                        for (auto [dx, dy] : deltas) {
                            int nx = x + dx, ny = y + dy;
                            if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
                                PType delta_p = old_p[x][y] - old_p[nx][ny];
                                PType force = delta_p;
                                VType &contr = velocity.get(nx, ny, -dx, -dy);
                                if (contr * rho[(int) field[nx][ny]] >= force) {
                                    contr -= force / rho[(int) field[nx][ny]];
                                    continue;
                                }
                                force -= contr * rho[(int) field[nx][ny]];
                                contr = 0;
                                velocity.add(x, y, dx, dy, force / rho[(int) field[x][y]]);
                                p[x][y] -= force / dirs[x][y];
                                total_delta_p -= force / dirs[x][y];
                            }
                        }
                        // cout << "deltas\n";
                    }
                }

                // Make flow from velocities
                velocity_flow = DynamicVectorField<VFType>(N, M);;
                bool prop = false;
                do {
                    UT += 2;
                    prop = 0;
                    for (size_t x = 0; x < N; ++x) {
                        for (size_t y = 0; y < M; ++y) {
                            if (field[x][y] != '#' && last_use[x][y] != UT) {
                                auto [t, local_prop, _] = propagate_flow(x, y, 1);
                                if (t > 0) {
                                    prop = 1;
                                }
                            }
                        }
                    }
                } while (prop);

                // Recalculate p with kinetic energy
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] == '#')
                            continue;
                        for (auto [dx, dy] : deltas) {
                            VType old_v = velocity.get(x, y, dx, dy);
                            VType new_v = velocity_flow.get(x, y, dx, dy);
                            if (old_v > 0) {
                                // cout << new_v.v << " " << old_v.v << endl;
                                assert(new_v <= old_v);
                                velocity.get(x, y, dx, dy) = new_v;
                                PType force = (old_v - new_v) * rho[(int) field[x][y]];
                                if (field[x][y] == '.')
                                    force *= 0.8;
                                if (field[x + dx][y + dy] == '#') {
                                    // cout << "# ";
                                    p[x][y] += force / dirs[x][y];
                                    total_delta_p += force / dirs[x][y];
                                } else {
                                    // cout << "!# ";
                                    p[x + dx][y + dy] += force / dirs[x + dx][y + dy];
                                    total_delta_p += force / dirs[x + dx][y + dy];
                                }
                            }
                            // cout << "recalc\n";
                        }
                    }
                }

                UT += 2;
                prop = false;
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] != '#' && last_use[x][y] != UT) {
                            if (VType::random01(rnd) < move_prob(x, y)) {
                                prop = true;
                                propagate_move(x, y, true);
                            } else {
                                propagate_stop(x, y, true);
                            }
                        }
                    }
                }

                // cout << "Tick " << i << ":\n";
                // for (size_t x = 0; x < N; ++x) {
                //     for (size_t y = 0; y < M; ++y) {
                //         cout << p[x][y] << " ";
                //     }
                //     cout << endl;
                // }
                // cout << endl;
                if (prop) {
                    cout << "Tick " << i << ":\n";
                    for (size_t x = 0; x < N; ++x) {
                        cout << field[x] << "\n";
                    }
                }
            }
        }

    };
}