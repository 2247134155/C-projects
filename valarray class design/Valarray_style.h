// Valarray.h

/* Put your solution in this file, we expect to be able to use
 * your epl::valarray class by simply saying #include "Valarray.h"
 *
 * We will #include "Vector.h" to get the epl::vector<T> class
 * before we #include "Valarray.h". You are encouraged to test
 * and develop your class using std::vector<T> as the base class
 * for your epl::valarray<T>
 * you are required to submit your project with epl::vector<T>
 * as the base class for your epl::valarray<T>
 */

#ifndef _Valarray_h
#define _Valarray_h

#include "Vector.h"
#include "vector"
#include <complex>
#include <cmath>
#include <string>

using epl::vector; // during development and testing
// using epl::vector; // after submission

//using namespace std;


template <typename T, typename R> class vect_wrap;

template<typename T1, typename T2 = vector<T1>>
using  valarray = vect_wrap<T1, T2>;



template <typename T> struct My_traits {
};

template <> struct My_traits<void> {
    enum {
        INT = 1,
        FLOAT = 2,
        DOUBLE = 3
    };
};

template <> struct My_traits<int> {
    enum {
        SRank = My_traits<void>::INT,
        CRank = false,
        VRank = false
    };
};

template <> struct My_traits<float> {
    enum {
        SRank = My_traits<void>::FLOAT,
        CRank = false,
        VRank = false
    };
};

template <> struct My_traits<double> {
    enum {
        SRank = My_traits<void>::DOUBLE,
        CRank = false,
        VRank = false
    };
};

template <typename T> struct My_traits<std::complex<T>> {
    enum {
        SRank = My_traits<T>::SRank,
        CRank = true,
        VRank = false
    };
};

template <typename T, typename R> struct My_traits<vect_wrap<T, R>> {
    enum {
        SRank = My_traits<T>::SRank,
        CRank = My_traits<T>::CRank,
        VRank = true
    };
};

template <uint32_t rank>
struct SType;

template <> struct SType<1> { using type = int;  };
template <> struct SType<2> { using type = float; };
template <> struct SType<3> { using type = double; };

template <typename T, bool complexity>
struct CType;

template <typename T> struct CType<T, true> {
    using type = std::complex<T>;
};

template <typename T> struct CType<T, false> {
    using type = T;
};

template <typename T, bool complexity>
using CompType = typename CType<T, complexity>::type;

template<typename T1, typename T2>
struct c_type {
    static const uint64_t t1_rank = My_traits<T1>::SRank;
    static const uint64_t t2_rank = My_traits<T2>::SRank;
    static const uint64_t max_rank = t1_rank > t2_rank ? t1_rank : t2_rank;
    using stype = typename SType<max_rank>::type;
    
    static const bool complexity = My_traits<T1>::CRank || My_traits<T2>::CRank;
    using type = CompType<stype, complexity>;
};

template <typename To, typename From, bool is_to_complex, bool is_from_complex>
struct convert;

template<typename To, typename From>
struct convert<To, From, true, true> {
    From val;
    convert(From _v) { val = _v; }
    To getValue() {
        return To{static_cast<typename SType<My_traits<To>::SRank>::type>(val.real()), static_cast<typename SType<My_traits<To>::SRank>::type>(val.imag())};
    }
};

template<typename To, typename From>
struct convert<To, From, true, false> {
    From val;
    convert(From _v) { val = _v; }
    To getValue() {
        return To{static_cast<typename SType<My_traits<To>::SRank>::type>(val)};
    }
};

template<typename To, typename From>
struct convert<To, From, false, false> {
    From val;
    convert(From _v) { val = _v; }
    To getValue() {
        return static_cast<To>(val);
    }
};

template<typename To, typename From, typename R=convert<To, From, My_traits<To>::CRank, My_traits<From>::CRank>>
struct promote : public R {
    promote(From _v) : R(_v) {}
};

template<typename T>
struct scalar {
    T val;
    scalar(T _v) : val(_v) {}
    
    T operator[](uint64_t index) const {
        return val;
    }
    
    uint64_t size() const {
        return -1;
    }
};

template<typename OperationProxy, typename T>
class My_Iter {
public:
    
    using iterator = My_Iter<OperationProxy, T>;
    
    OperationProxy* proxy;
    uint64_t index;
    
    My_Iter(OperationProxy* _p, int _i) : proxy(_p), index(_i) {}
    
    T operator*() {
        return proxy->operator[](index);
    }
    
    
    My_Iter& operator++() {
        ++index;
        return *this;
    }
    
    
    My_Iter operator++(int) {
        iterator t{*this};
        ++index;
        return t;
    }
    
    
    My_Iter operator+(int64_t n) {
        return iterator{proxy, index + n};
    }
    
   
    My_Iter operator-(int64_t n) {
        return iterator{proxy, index - n};
    }
    
    
    My_Iter& operator+=(int64_t n) {
        index += n;
        return *this;
    }
    
    
    My_Iter& operator-=(int64_t n) {
        index -= n;
        return *this;
    }
    
    
    T operator[](int64_t n) {
        return proxy->operator[](index + n);
    }
    
    bool operator == (const My_Iter& rhs) const {
        return this->index == rhs.index ? true : false;
    }
    
};


template <typename L, typename R, typename FstOper>
struct OneProxy {
    using RetType = typename FstOper::result_type;
    
    const vect_wrap<L, R>& val;
    FstOper op;
    
    OneProxy(const vect_wrap<L, R>& _v, FstOper _op) : val(_v), op(_op) {}
    
    RetType operator[](uint64_t index) const {
        return op(this->val[index]);
    }
    
    uint64_t size() const { return (uint64_t)val.size(); }
    
    using iterator = My_Iter<OneProxy, L>;
    using const_iterator = My_Iter<OneProxy, const L>;
    
    iterator begin() {
        return iterator(this, 0);
    }
    
    const_iterator begin() const {
        return const_iterator(this, 0);
    }
    
    iterator end() {
        return iterator(this, this->size());
    }
    
    const_iterator end() const {
        return const_iterator(this, this->size());
    }
};

template <typename Operator, typename LRetType, typename LType, typename RRetType, typename RType>
struct PareProxy {
    using RetType = typename c_type<LRetType, RRetType>::type;
    
    const vect_wrap<LRetType, LType>& left;
    const vect_wrap<RRetType, RType>& right;
    Operator op;
    
    
    PareProxy(Operator _op, const vect_wrap<LRetType, LType>& _l, const vect_wrap<RRetType, RType>& _r) : op(_op), left(_l), right(_r) {}
    
    RetType operator[](uint64_t index) const {
        auto lValue = promote<RetType, LRetType>(this->left[index]);
        auto rValue = promote<RetType, RRetType>(this->right[index]);
        
        return this->op(lValue.getValue(), rValue.getValue());
    }
    
    
    uint64_t size() const {
        return std::min((uint64_t)left.size(), (uint64_t)right.size());
    }
};


template <typename Operator, typename LRetType, typename LType, typename RRetType>
struct PareProxy<Operator, LRetType, LType, RRetType, scalar<RRetType>> {
    using RetType = typename c_type<LRetType, RRetType>::type;
    
    
    const vect_wrap<LRetType, LType>& left;
    scalar<RRetType> right;
    Operator op;
    
    
    bool order;
    
    PareProxy(Operator _op, const vect_wrap<LRetType, LType>& _l, scalar<RRetType> _r, bool _b) : op(_op), left(_l), right(_r), order(_b) {}
    
    RetType operator[](uint64_t index) const {
        auto lValue = promote<RetType, LRetType>(this->left[index]);
        auto rValue = promote<RetType, RRetType>(this->right[index]);
        
        if (order) {
            return this->op(lValue.getValue(), rValue.getValue());
        } else {
            return this->op(rValue.getValue(), lValue.getValue());
        }
    }
    
    
    uint64_t size() const {
        return (uint64_t)left.size();
    }
};


template <typename Operation, typename LRetType, typename LType, typename RRetType=LRetType, typename RType=scalar<RRetType>, typename R=PareProxy<Operation, LRetType, LType, RRetType, RType>>
struct Proxy : public R {
    using RetType = typename c_type<LRetType, RRetType>::type;
    using iterator = My_Iter<Proxy, RetType>;
    using const_iterator = My_Iter<Proxy, const RetType>;
    Proxy(Operation _op, const vect_wrap<LRetType, LType>& _l, const vect_wrap<RRetType, RType>& _r) : R(_op, _l, _r) {}
    Proxy(Operation _op, const vect_wrap<LRetType, LType>& _l, scalar<RRetType> _r, bool order) : R(_op, _l, _r, order) {}
    
    iterator begin() {
        return iterator(this, 0);
    }
    
    const_iterator begin() const {
        return const_iterator(this, 0);
    }
    
    iterator end() {
        return iterator(this, this->size());
    }
    
    const_iterator end() const {
        return const_iterator(this, this->size());
    }
    
};

template <typename Arg, typename Result>
struct square_root : public std::unary_function<Arg, Result>  {
    Result operator() (const Arg& arg) const {
        return std::sqrt(arg);
    }
};

template <typename T, typename R=vector<T> >
struct vect_wrap : public R {
public:
    static const int zyh = 1;
    vect_wrap() : R() {}
    
    explicit vect_wrap(uint64_t size) : R(size) {}
    
    vect_wrap(std::initializer_list<T> il) : R(il) {}
    
    template <typename Operator, typename T1, typename R1, typename T2, typename R2>
    vect_wrap(Operator op, const vect_wrap<T1, R1>& x, const vect_wrap<T2, R2>& y) : R(op, x, y) {}
    
    template <typename Operator, typename T1, typename R1, typename T2>
    vect_wrap(Operator op, const vect_wrap<T1, R1>& x, scalar<T2> y, bool order=true) : R(op, x, y, order) {}
    
    
    template <typename L1, typename R1, typename FstOper>
    vect_wrap(const vect_wrap<L1, R1>& rhs, FstOper op) : R(rhs, op) {}
    
    
    template <typename L1, typename R1>
    vect_wrap(const vect_wrap<L1, R1>& rhs) {
        for (int i=0; i < rhs.size(); ++i) {
            this->push_back((T)rhs[i]);
        }
    }
    
    vect_wrap<T, R>& operator=(const vect_wrap<T, R>& that) {
        std::cout<<"11111111"<<std::endl;
        if (this != &that) {
            int total_size = this->size();
            int temp = 0;
            while(temp != total_size) {
                ++temp;
                this->pop_back();
            }
            for (int i=0; i < that.size(); ++i) {
                this->push_back(that[i]);
            }
        }
        return *this;
    }
    
    
    template <typename L1, typename R1>
    vect_wrap<T, R>& operator=(const vect_wrap<L1, R1>& that) {
        
        
        int total_size = this->size();
        int record_size = that.size();
        int temp = 0;
        
        for(unsigned int i = 0; i < that.size(); i++){
            (*this)[i] = (T)that[i];
        }
        return *this;
        //return *this;
    }
    
    
    vect_wrap<T, R>& operator=(const T& val) {
        for (int i=0; i < this->size(); ++i) {
            this->operator[](i) = val;
        }
    }
    
    template<class Operator>
    T accumulate(Operator op, T init) {
        for (int i=0; i < this->size();  ++i) {
            init = op(init, this->operator[](i));
        }
        return init;
    }
    
    template<class Operator>
    T accumulate(Operator op) {
        T init = this->operator[](0);
        for (int i=1; i < this->size();  ++i) {
            init = op(init, this->operator[](i));
        }
        return init;
    }
    T sum() {
        T result{};
        return accumulate(std::plus<T>(), result);
    }
    
    template<typename FstOper>
    vect_wrap<typename FstOper::result_type, OneProxy<T, R, FstOper>>
    apply(FstOper op) const {
        return vect_wrap<typename FstOper::result_type, OneProxy<T, R, FstOper>>(*this, op);
    }
    
    template <typename T1=typename c_type<T, CompType<double, My_traits<T>::CRank> >::type>
    vect_wrap<T1, OneProxy<T, R, square_root<T, T1>>>
    sqrt() {
        return apply(square_root<T, T1>());
    }
};


template <typename T1, typename R1>
vect_wrap<T1, OneProxy<T1, R1, std::negate<T1>> >
operator-(const vect_wrap<T1, R1>& x) {
    return x.apply(std::negate<T1>());
}


template <template<typename> class OP>
struct MyStruct{
    
    template <typename T1,typename R1 ,typename T2 ,typename R2,typename RetType=typename c_type<T1, T2>::type>
    static vect_wrap<RetType,Proxy<OP<RetType>,T1,R1,T2,R2>> foo(const vect_wrap<T1, R1>& x, const vect_wrap<T2, R2>& y){
        return vect_wrap<RetType, Proxy<OP<RetType>, T1, R1, T2, R2>>(OP<RetType>(), x, y);
    }
    
    template <typename T1, typename R1, typename T2, typename RetType=typename c_type<T1, T2>::type>
    static vect_wrap<RetType, Proxy<OP<RetType>, T1, R1, T2> > foo(const vect_wrap<T1, R1>& x, T2 y){
        return vect_wrap<RetType, Proxy<OP<RetType>, T1, R1, T2>>(OP<RetType>(), x, scalar<T2>(y), true);
    }
    
    template <typename T1, typename T2, typename R2, typename RetType=typename c_type<T1, T2>::type>
    static vect_wrap<RetType, Proxy<OP<RetType>, T2, R2, T1> > foo(T1 y, const vect_wrap<T2, R2>& x){
        return vect_wrap<RetType, Proxy<OP<RetType>, T2, R2, T1>>(OP<RetType>(), x, scalar<T1>(y), false);
    }
    
};
template <typename T1 ,typename T2>
auto operator +(const T1& x , const T2& y)->decltype(MyStruct<std::plus>::foo(x,y)){
    return MyStruct<std::plus>::foo(x,y);
}
template <typename T1 ,typename T2>
auto operator -(const T1& x , const T2& y)->decltype(MyStruct<std::minus>::foo(x,y)){
    return MyStruct<std::minus>::foo(x,y);
}
template <typename T1 ,typename T2>
auto operator *(const T1& x , const T2& y)->decltype(MyStruct<std::multiplies>::foo(x,y)){
    return MyStruct<std::multiplies>::foo(x,y);
}

template <typename T1 ,typename T2>
auto operator /(const T1& x , const T2& y)->decltype(MyStruct<std::divides>::foo(x,y)){
    return MyStruct<std::divides>::foo(x,y);
}


/////////end





//template <typename T1, typename R1, typename T2, typename R2, typename RetType=typename c_type<T1, T2>::type>
//vect_wrap<RetType, Proxy<std::plus<RetType>, T1, R1, T2, R2> >
//operator+(const vect_wrap<T1, R1>& x, const vect_wrap<T2, R2>& y) {
//
//    return vect_wrap<RetType, Proxy<std::plus<RetType>, T1, R1, T2, R2>>(std::plus<RetType>(), x, y);
//}
//
//template <typename T1, typename R1, typename T2, typename RetType=typename c_type<T1, T2>::type>
//vect_wrap<RetType, Proxy<std::plus<RetType>, T1, R1, T2> >
//operator+(const vect_wrap<T1, R1>& x, T2 y) {
//    //cout<<"2+"<<endl;
//    return vect_wrap<RetType, Proxy<std::plus<RetType>, T1, R1, T2>>(std::plus<RetType>(), x, scalar<T2>(y), true);
//}
//
//template <typename T1, typename T2, typename R2, typename RetType=typename c_type<T1, T2>::type>
//vect_wrap<RetType, Proxy<std::plus<RetType>, T2, R2, T1> >
//operator+(T1 y, const vect_wrap<T2, R2>& x) {
//    //cout<<"3+"<<endl;
//    return vect_wrap<RetType, Proxy<std::plus<RetType>, T2, R2, T1>>(std::plus<RetType>(), x, scalar<T1>(y), false);
//}
//
//template <typename T1, typename R1, typename T2, typename R2, typename RetType=typename c_type<T1, T2>::type>
//vect_wrap<RetType, Proxy<std::minus<RetType>, T1, R1, T2, R2> >
//operator-(const vect_wrap<T1, R1>& x, const vect_wrap<T2, R2>& y) {
//    //cout<<"1-"<<endl;
//    return vect_wrap<RetType, Proxy<std::minus<RetType>, T1, R1, T2, R2>>(std::minus<RetType>(), x, y);
//}
//
//template <typename T1, typename R1, typename T2, typename RetType=typename c_type<T1, T2>::type>
//vect_wrap<RetType, Proxy<std::minus<RetType>, T1, R1, T2>>
//operator-(const vect_wrap<T1, R1>& x, T2 y) {
//    //cout<<"2-"<<endl;
//    return vect_wrap<RetType, Proxy<std::minus<RetType>, T1, R1, T2>>(std::minus<RetType>(), x, scalar<T2>(y), true);
//}
//
//template <typename T1, typename T2, typename R2, typename RetType=typename c_type<T1, T2>::type>
//vect_wrap<RetType, Proxy<std::minus<RetType>, T2, R2, T1> >
//operator-(T1 y, const vect_wrap<T2, R2>& x) {
//    //cout<<"3-"<<endl;
//    return vect_wrap<RetType, Proxy<std::minus<RetType>, T2, R2, T1>>(std::minus<RetType>(), x, scalar<T1>(y), false);
//}
//
//template <typename T1, typename R1, typename T2, typename R2, typename RetType=typename c_type<T1, T2>::type>
//vect_wrap<RetType, Proxy<std::multiplies<RetType>, T1, R1, T2, R2> >
//operator*(const vect_wrap<T1, R1>& x, const vect_wrap<T2, R2>& y) {
//    //cout<<"1*"<<endl;
//    return vect_wrap<RetType, Proxy<std::multiplies<RetType>, T1, R1, T2, R2>>(std::multiplies<RetType>(), x, y);
//}
//
//template <typename T1, typename R1, typename T2, typename RetType=typename c_type<T1, T2>::type>
//vect_wrap<RetType, Proxy<std::multiplies<RetType>, T1, R1, T2>>
//operator*(const vect_wrap<T1, R1>& x, T2 y) {
//    //cout<<"2*"<<endl;
//    return vect_wrap<RetType, Proxy<std::multiplies<RetType>, T1, R1, T2>>(std::multiplies<RetType>(), x, scalar<T2>(y), true);
//}
//
//template <typename T1, typename T2, typename R2, typename RetType=typename c_type<T1, T2>::type>
//vect_wrap<RetType, Proxy<std::multiplies<RetType>, T2, R2, T1> >
//operator*(T1 y, const vect_wrap<T2, R2>& x) {
//    //cout<<"3*"<<endl;
//    return vect_wrap<RetType, Proxy<std::multiplies<RetType>, T2, R2, T1>>(std::multiplies<RetType>(), x, scalar<T1>(y), false);
//}
//
//template <typename T1, typename R1, typename T2, typename R2, typename RetType=typename c_type<T1, T2>::type>
//vect_wrap<RetType, Proxy<std::divides<RetType>, T1, R1, T2, R2> >
//operator/(const vect_wrap<T1, R1>& x, const vect_wrap<T2, R2>& y) {
//   
//    return vect_wrap<RetType, Proxy<std::divides<RetType>, T1, R1, T2, R2>>(std::divides<RetType>(), x, y);
//}
//
//template <typename T1, typename R1, typename T2, typename RetType=typename c_type<T1, T2>::type>
//vect_wrap<RetType, Proxy<std::divides<RetType>, T1, R1, T2>>
//operator/(const vect_wrap<T1, R1>& x, T2 y) {
//    return vect_wrap<RetType, Proxy<std::divides<RetType>, T1, R1, T2>>(std::divides<RetType>(), x, scalar<T2>(y), true);
//}
//
//template <typename T1, typename T2, typename R2, typename RetType=typename c_type<T1, T2>::type>
//vect_wrap<RetType, Proxy<std::divides<RetType>, T2, R2, T1> >
//operator/(T1 y, const vect_wrap<T2, R2>& x) {
//    return vect_wrap<RetType, Proxy<std::divides<RetType>, T2, R2, T1>>(std::divides<RetType>(), x, scalar<T1>(y), false);
//}

template<typename T, typename R>
std::ostream& operator<<(std::ostream& out, const vect_wrap<T, R>& val) {
    const char* p = "";
    for (int i=0; i < val.size(); ++i) {
        out << p << val.operator[](i);
        p = ", ";
    }
    out << std::endl;
    return out;
}
#endif /* _Valarray_h */