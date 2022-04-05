#ifndef MATRIX_SRC_MATRIX_H
#define MATRIX_SRC_MATRIX_H

#include <algorithm>
#include <iostream>
template <typename T, size_t H, size_t W>
class matrix {
private:
    T** array;
    using matrix_minor = matrix<T, std::max<size_t>(H - 1, 1), std::max<size_t>(W - 1, 1)>;

public:

    explicit matrix(): array(new T* [H]) {
        for ( size_t i = 0; i < H; i++){
            array[i] = new T[W];
        }
        for( size_t i = 0;i < H;i++)
            for(size_t j=0; j < W;j++)
                array[i][j]=0;
    }
    ~matrix() ;
    matrix (const T& value): array( new T* [H]){
        for ( size_t i = 0; i < H; i++){
            array[i] = new T[W];
        }
        for(size_t i=0; i<H;i++)
            for(size_t j=0; j<W;j++)
                array[i][j]=value;
    }
    matrix(const matrix& other): array(new T* [H]){
        for ( size_t i = 0; i < H; i++){
            array[i] = new T[W];
        }
        for(size_t i=0; i < H;i++)
            for(size_t j=0;j<W;j++)
                array[i][j]=other.array[i][j];
    }

    const T& at(const size_t& i, const size_t& j) const{
        return array[i][j];
    }

    T& at(const size_t& i, const size_t& j){
        return array[i][j];
    }
    matrix& operator = (const matrix& other){
        if(this!=&other){
            for(size_t i=0;i<H;i++)
                for(size_t j=0;j<W;j++)
                    array[i][j]=other.array[i][j];
        }
        return *this;
    }
    const matrix operator+ () const {
        return *this;
    }
    matrix& operator *=(const matrix<T,W,W>& other){
        matrix temp = (*this);
        for(size_t i=0;i<H;i++ ) {
            for (size_t j = 0; j < W; j++) {
                array[i][j] = 0;
            }
        }

        for(size_t i=0;i<H;i++) {
            for (size_t j = 0; j < W; j++) {
                for (size_t l = 0; l < W; l++){
                    array[i][j] += temp.array[i][l] * other.at(l,j);}
            }

        }
        return *this;

    }


    const matrix operator- () const {
        matrix temp(*this);
        for(size_t i=0;i<H;i++)
            for(size_t j=0;j<W;j++)
                temp.array[i][j]=-array[i][j];

        return temp;
    }
    bool  operator != (  const matrix& other ) const {
        for ( size_t i = 0; i < H; ++i){
            for ( size_t j = 0; j < W; ++j){
                if ( array[i][j] != other.array[i][j] ){
                    return true;
                }
            }
        }
        return false;
    }
    bool operator == ( const matrix& other) const {
        for ( size_t i = 0; i < H; ++i){
            for ( size_t j = 0; j < W; ++j){
                if ( array[i][j] != other.array[i][j] ){
                    return false;
                }
            }
        }
        return true;
    }
    matrix& operator *= (const T& value){
        for(size_t i=0;i<H;i++)
            for(size_t j=0;j<W;j++)
                array[i][j]*=value;

        return *this;
    }
    matrix& operator += ( const matrix& other )  {
        for (size_t i = 0; i < H; i++) {
            for (size_t j = 0; j < W; j++) {
                array[i][j] += other.array[i][j];
            }
        }
        return *this;
    }

    matrix& operator -= ( const matrix& other)  {
        for (size_t i = 0; i < H; i++) {
            for (size_t j = 0; j < W; j++) {
                array[i][j] -= other.array[i][j];
            }
        }
        return *this;
    }

    matrix operator + ( const matrix& other) const{
        matrix result = *this;
        result += other;
        return result;
    }
    matrix operator - ( const matrix& other) const{
        matrix result = *this;
        result -= other;
        return result;
    }
    matrix<T,W,H> transposed() const {
        matrix<T,W,H> temp;
        for(size_t i=0;i<H;i++)
            for(size_t j=0;j<W;j++)
                temp.at(j,i) = array[i][j];
        return temp;
    }
    T trace() const{
        T trace_ = T();
        for ( size_t i = 0; i < W; i++){
            trace_ += array[i][i];
        }
        return trace_;
    }
    T det() const{
        T res=0;
        if( W==1 )
            return array[0][0];

        for(size_t k=0; k<W;k++) {
            matrix_minor matrix_(0);

            for(size_t i=0;i<W-1;i++)
                for(size_t j=0;j<k;j++)
                    matrix_.at(i,j)=array[i+1][j];

            for(size_t i=0;i<W-1;i++)
                for(size_t j=k;j<W+1;j++)
                    matrix_.at(i,j)=array[i+1][j+1];


            if ( k % 2 == 0){
                res += array[0][k] * matrix_.det();
            }
            else {
                res += (-1)*array[0][k] * matrix_.det();
            }
        }
        return res;

    }

};

template<typename T,size_t H, size_t W>
matrix<T,H,W> operator+( const matrix<T,H,W> first, const T value ){
    return matrix<T,H,W>(first) += value;
}

template<typename T,size_t H, size_t W>
matrix<T,H,W> operator-(const matrix<T,H,W>& first, const T& value){
    return matrix<T,H,W>(first) -= value;
}

template<typename T,size_t H, size_t W>
matrix<T,H,W> operator+(const T& value,const matrix<T,H,W>& first) {
    return matrix<T,H,W>(first) += value;
}
template<typename T,size_t H, size_t W>
matrix<T,H,W> operator-(const T& value,const matrix<T,H,W>& first){
    return matrix<T,H,W>(first) -= value;
}
template<typename T,size_t H, size_t W>
matrix<T,H,W> operator*(const T& value,const matrix<T,H,W>& first){
    return matrix<T,H,W>(first)*=value;
}

template<typename T,size_t H, size_t W>
matrix<T,H,W> operator*(const matrix<T,H,W>& first, const T& value) {
    return matrix<T, H, W>(first) *= value;
}
template<typename T,size_t H, size_t W,size_t K>
matrix<T,H,K> operator*(const matrix<T,H,W>& first, const matrix<T,W,K>& second) {
    return matrix<T,H,W>(first) *= second;

}
#endif /// MATRIX_SRC_MATRIX_H.
