#include <bits/stdc++.h>

using namespace std;

template<int32_t N, int32_t K, bool fast = false>
requires(N <= 8 * sizeof(intmax_t) && N > 0 && K > 0)
struct Fixed
{
    using type = std::conditional_t<
        fast == false,
        std::conditional_t<N <= 8, int8_t,
            std::conditional_t<N <= 16, int16_t,
                std::conditional_t<N <= 32, int32_t,
                    std::conditional_t<N <= 64, int64_t, void>>>>,
        std::conditional_t<N <= 8, int_fast8_t,
                std::conditional_t<N <= 16, int_fast16_t,
                        std::conditional_t<N <= 32, int_fast32_t,
                                std::conditional_t<N <= 64, int_fast64_t, void>>>>>;
    type v;

    static Fixed random01(auto rnd)
    {
        return Fixed<N, K, fast>::from_raw((rnd() & ((1 <<  K) - 1)));
    }

    constexpr Fixed() : v(0)
    {
        assert(N <= 8 * sizeof(intmax_t));
    }
    constexpr Fixed(int v_) : Fixed()
    {
        v = v_ << K;
    }
    constexpr Fixed(float f): Fixed()
    {
        v = f * (1 << K);
    }
    constexpr Fixed(double f): Fixed()
    {
        v = f * (1 << K);
    }


    template<int N_, int K_, bool fast_ = false>
    constexpr Fixed(const Fixed<N_, K_, fast_>& fx)
    {
        if constexpr (K_ > K)
        {
            v = (intmax_t)fx.v >> (K_ - K);
        }
        else
        {
            v = (intmax_t)fx.v << (K - K_);
        }
    }

    template<int N_, int K_, bool fast_ = false>
    explicit operator Fixed<N_, K_, fast_>()
    {
        return Fixed<N_, K_, fast_>(*this);
    }

    static constexpr Fixed from_raw(intmax_t x) {
        Fixed ret{};
        ret.v = x;
        return ret;
    }

    const Fixed &add_impl(const Fixed &b) {
        v += b.v;
        return *this;
    }
    template<int N_, int K_, bool fast_ = false>
    friend Fixed operator+(Fixed a, Fixed<N_, K_, fast_> b) {
        return a.add_impl(b);
    }
    friend Fixed operator+(Fixed a, Fixed b) {
        return Fixed::from_raw(a.v + b.v);
    }

    const Fixed &sub_impl(const Fixed &b) {
        v -= b.v;
        return *this;
    }
    template<int N_, int K_, bool fast_ = false>
    friend Fixed operator-(Fixed a, Fixed<N_, K_, fast_> b) {
        return a.sub_impl(b);
    }
    friend Fixed operator-(Fixed a, Fixed b) {
        return Fixed::from_raw(a.v - b.v);
    }

    const Fixed &mult_impl(const Fixed &b) {
        v *= b.v;
        v = v >> K;
        return *this;
    }
    template<int N_, int K_, bool fast_ = false>
    friend Fixed operator*(Fixed a, Fixed<N_, K_, fast_> b) {
        return a.mult_impl(b);
    }
    friend Fixed operator*(Fixed a, Fixed b) {
        return Fixed::from_raw(((intmax_t) a.v * b.v) >> K);
    }

    const Fixed &div_impl(const Fixed &b) {
        v = v << K;
        v /= b.v;
        return *this;
    }
    template<int N_, int K_, bool fast_ = false>
    friend Fixed operator/(Fixed a, Fixed<N_, K_, fast_> b) {
        return a.div_impl(b);
    }
    friend Fixed operator/(Fixed a, Fixed b) {
        return Fixed::from_raw(((intmax_t) a.v << K) / b.v);
    }

    friend Fixed &operator+=(Fixed &a, Fixed b) {
        return a = a + b;
    }

    friend Fixed &operator-=(Fixed &a, Fixed b) {
        return a = a - b;
    }

    friend Fixed &operator*=(Fixed &a, Fixed b) {
        return a = a * b;
    }

    friend Fixed &operator/=(Fixed &a, Fixed b) {
        return a = a / b;
    }

    friend Fixed operator-(Fixed x) {
        return Fixed::from_raw(-x.v);
    }

    Fixed abs(Fixed x) {
        if (x.v < 0) {
            x.v = -x.v;
        }
        return x;
    }

    friend ostream &operator<<(ostream &out, Fixed x) {
        return out << (double)(x.v / (double) ((intmax_t)1 << K));
    }

    explicit constexpr operator float() { return float(v) / ((intmax_t)1 << K); }
    explicit constexpr operator double() { return double(v) / ((intmax_t)1 << K); }
    auto operator<=>(const Fixed&) const = default;
    bool operator==(const Fixed&) const = default;

    bool more_impl(const Fixed &b) {
        return v > b.v;
    }
    template<int N_, int K_, bool fast_ = false>
    friend bool operator>(Fixed a, Fixed<N_, K_, fast_> b) {
        return a.more_impl(b);
    }
    bool operator>(const Fixed& cmp) const {
        return v > cmp.v;
    }

    bool less_impl(const Fixed &b) {
        return v < b.v;
    }
    template<int N_, int K_, bool fast_ = false>
    friend bool operator<(Fixed a, Fixed<N_, K_, fast_> b) {
        return a.less_impl(b);
    }
    bool operator<(const Fixed& cmp) const {
        return v < cmp.v;
    }

    bool moreeq_impl(const Fixed &b) {
        return v >= b.v;
    }
    template<int N_, int K_, bool fast_ = false>
    friend bool operator>=(Fixed a, Fixed<N_, K_, fast_> b) {
        return a.moreeq_impl(b);
    }
    bool operator>=(const Fixed& cmp) const {
        return v >= cmp.v;
    }

    bool lesseq_impl(const Fixed &b) {
        return v < b.v;
    }
    template<int N_, int K_, bool fast_ = false>
    friend bool operator<=(Fixed a, Fixed<N_, K_, fast_> b) {
        return a.lesseq_impl(b);
    }
    bool operator<=(const Fixed& cmp) const {
        return v <= cmp.v;
    }
};