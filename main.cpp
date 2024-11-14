#include <iostream>
#include <string>
#include <memory>
#include <math.h>


enum DirectionCol {Ia, NoIa}; 
enum DirectionRow {Rerun, Endrun};
enum Stage {FIRSTSTAGE, OTHERSTAGE};
void firstPhase(const std::shared_ptr<int[]> &R,
        const std::shared_ptr<int[]>& Q,
        const std::shared_ptr<int[]>& u,
        const std::shared_ptr<int[]>& v,
        const std::shared_ptr<int[]>& ik,
        const std::shared_ptr<int[]>& jk,
        const std::shared_ptr<int[]>& epsilonRow,
        int n,
        Stage stage);
void fillCustom(const std::shared_ptr<int[]> & arr, int size, int value){
    for(int i = 0; i < size; i++) {
        arr[i] = value; 
    }
}
void fillCustom(const std::shared_ptr<int[]> & arr, int size, const std::shared_ptr<int[]> &src){
    for(int i = 0; i < size; i++) {
        arr[i] = src[i]; 
    }
}
//R, u, v, i, j, size, true
int q(const std::shared_ptr<int[]>& R,  const std::shared_ptr<int[]> &u, const std::shared_ptr<int[]> &v, int i, int j,  int size, bool reversed){
    int n = sqrt(size);
    if(reversed){
        int tmp = i; 
        i = j; 
        j = tmp; 
    }
    return ((u[i] + v[j]) ==  R[ i * n + j ]) ? 1 : 0; 
}
void printSquaredMatrix(std::string name, const std::shared_ptr<int[]> & arr, int size){
    std::cout << "********Matrix display of " << name << " : ***** " << std::endl;
    int n = sqrt(size); 
    for(int i = 0; i < size; i++){
        std::cout << " " << arr[i] << " "; 
        if(i%n == n-1){
            std::cout << "" << std::endl ;
        }
    }
}


void initialization(const std::shared_ptr<int[]> &R, const std::shared_ptr<int[]> &Q, const std::shared_ptr<int[]> &u, const std::shared_ptr<int[]> &v, int size, Stage stage){
    // ai  = max j rij , for each i 
    // so we traverse the matrix line by line 
    // and for each line we check which ele has the biggest weight
    int n = sqrt(size);
    bool inversed = false;
    if(stage == FIRSTSTAGE){

        // to calculate the sum of max for the columns.
        const std::shared_ptr<int[]> a(new int[n]);
        const std::shared_ptr<int[]> b(new int[n]);

        int aS = 0, ai = 0; 
        int bS = 0;
        for(int i = 0; i < size ; i++){
            if(i/n == 0){
                v[i] = R[i];  
                u[i/n] = R[i];
            }
            if(R[i] > u[i/n]) {
                u[i/n] = R[i]; 
            }
            if(R[i] > v[i%n]){
                v[i%n] = R[i];  
            }
            if(i%n == n-1){
                aS+= u[i/n]; 
            }
            if(i/n == n-1){
                bS += v[i%n];             
            }
        }
        if (aS > bS) {
            fillCustom(u, n, 0);
            inversed = true;
        }else {
            fillCustom(v, n, 0);  
            inversed = false; 
        }

    } else {
        // TODO this might be a mistake;
        inversed = false;  
    }
    bool first = true; 
    for(int i = 0 ; i < size ; i++){
        Q[i] = q(R, u, v, i/n, i%n, size, inversed);  
        if(Q[i] == 1 && first){
            bool falseAlarm = false; 
            int j = i%n; 
            for (int l = 0; l < n ; l++) {
                if(Q[j + l*n] == 2) {
                    falseAlarm = true;   
                }
            }
            if(!falseAlarm){
                Q[i] = 2;  
                first = false;
            }
        }
        if(i%n == n-1) {
            first = true; 
        }
    }
    printSquaredMatrix("Q in init ", Q, size);
    return;
}
/**
 * ik: pointer to list of i indices with k values as key.
 * jk: pointer to the list of j indices with k values as key.
 * Q: pointer to the associated Q matrix of the paper.  
 * i: index of the row i.
 * n : the size of the row/col
 */
DirectionCol transferCol(
        const std::shared_ptr<int[]>& ik,
        const std::shared_ptr<int[]>& jk,
        const std::shared_ptr<int[]> &Q,
        const std::shared_ptr<int[]> &epsilon,
        const std::shared_ptr<int> &k,
        int i, 
        int n){

    for(int h = 0 ; h < n ; h++){
        if(Q[i*n + h] == 2) {
            return NoIa;
        }
    }
    Q[ik[*k]*n + jk[*k - 1]] = 2;
    ik[*k] = 0; 
    // reverse the alternative path ( is it a augmented path ? ), so we are gaining one edge .
    while(*k > 1){
        Q[ik[*k-1]*n + jk[*k - 1]] = 1;
        (*k)--; 
        Q[ik[*k]*n + jk[*k - 1]] = 2;
        ik[*k] = 0; 
    }
    epsilon[*k-1] = 0; 
    while(*k < n){
        (*k)++; 
        // TODO i am not sure if it is a mistake from the paper or if i should do it an other way 
        epsilon[*k-1] = 0;     
    }
    *k = 1;
    // clean up part in case we had constructed some essential rows. 
    // I don't see it's necessity in the implementation, as i assume 
    // when we call again the first stage, it is reset. 
    return Ia; 
}

/**
 * ik: pointer to list of i indices with k values as key.
 * jk: pointer to the list of j indices with k values as key.
 * Q: pointer to the associated Q matrix of the paper.  
 * epsilon: pointer to the record of essential rows.  
 * k: pointer to the int representing Tally of length of seq of 1's and 1*'s. 
 * j: index of the col j.
 * n : the size of the row/col.
 */
DirectionRow transferRow(
        const std::shared_ptr<int[]>& ik,
        const std::shared_ptr<int[]>& jk,
        const std::shared_ptr<int[]>& Q,
        const std::shared_ptr<int[]>& epsilon,
        const std::shared_ptr<int>& k,
        int j,
        int n){
    int i = 0 ; 
    bool firstPass = true; 
    do {
        if(!firstPass) {
            (*k)--; 
        } else {
            firstPass = false; 
        }
        while(i < n){
            if(Q[i*n + j] == 1){
                (*k)++; 
                for(int l = 0; l < n; l++){
                    if(i == ik[l]) {
                        break;
                    }
                }
                return Rerun; 
            }
            i++; 
        }

        if(epsilon[ik[*k]] <=0){
            epsilon[ik[*k]] = 1; 
            i = ik[*k];
            j = jk[*k-1];
            ik[*k] = 0; 
            jk[*k-1] = 0; 
        }
    }while(*k > 1); 

    return Endrun;



}

/*
 * Q: pointer to the associated Q matrix of the paper.  
 * epsilon: pointer to the record of essential rows.  
 * i: index of the row i.
 * j: index of the col j.
 * k: pointer to the int representing Tally of length of seq of 1's and 1*'s. 
 * n: size of rows/cols. (nxn) matrix.
 */
void transfer(
        const std::shared_ptr<int[]> &R,
        const std::shared_ptr<int[]>& Q,
        const std::shared_ptr<int[]>& epsilon,
        const std::shared_ptr<int[]>& u,
        const std::shared_ptr<int[]>& v,
        const std::shared_ptr<int[]>& ik,
        const std::shared_ptr<int[]>& jk,
        int i,
        int j,
        const std::shared_ptr<int> &k,
        int n){
    // shift chunk 
    // the column part
    ik[1] = i; 
    jk[0] = j;
    //const std::shared_ptr<int[]>& ik,
    //const std::shared_ptr<int[]>& jk,
    //const std::shared_ptr<int[]> &Q,
    //const std::shared_ptr<int[]> &epsilon,
    //int* k,
    //int i, 
    //int n){
    DirectionCol d = transferCol(ik, jk, Q, epsilon, k, i, n); 
    if(d == NoIa){
        transferRow(ik, jk,  Q, epsilon, k, j, n); 
        return; 
    }else {
        firstPhase(R, Q, u, v, ik, jk, epsilon,  n, OTHERSTAGE);
        return;  
    }
}

/**
 * Q: pointer to the associated Q matrix of the paper.  
 * n: size of rows/cols. (nxn) matrix.
 */
void firstPhase(
        const std::shared_ptr<int[]> &R,
        const std::shared_ptr<int[]>& Q,
        const std::shared_ptr<int[]>& u,
        const std::shared_ptr<int[]>& v,
        const std::shared_ptr<int[]>& ik,
        const std::shared_ptr<int[]>& jk,
        const std::shared_ptr<int[]>& epsilonRow,
        int n,
        Stage stage){
    initialization(R, Q, u, v, n*n, stage);
    const std::shared_ptr<int> k = std::make_unique<int>(); 
    *k = 1;
    for(int j=0; j < n; j++){
        for(int i=0; i < n; i++){
            if(Q[i*n + j] == 2){
                // break the row loop.
                break; 
            }
        }
        for(int i=0; i < n; i++){
            if(Q[i*n + j] == 1){
                transfer(R, Q, epsilonRow, u, v, ik, jk, i, j, k, n);
            }
        }
    }
    printSquaredMatrix("Q after first Phase ", Q, n*n);
}

/**
 * R: pointer to the input matrix R(nxn) 
 * u: pointer to an array representing the cover of rows.
 * v: pointer to an array representing the cover of columns. 
 * epsilonRow: pointer to the record of essential rows.  
 * epsilonCol: pointer to the record of essential cols.  
 * n: size of rows/cols. (nxn) matrix.
 */
void secondPhase(const std::shared_ptr<int[]>& R, const std::shared_ptr<int[]>& Q,  const std::shared_ptr<int[]>& u, const std::shared_ptr<int[]>& v, const std::shared_ptr<int[]>& ik, const std::shared_ptr<int[]>& jk, const std::shared_ptr<int[]>& epsilonRow, const std::shared_ptr<int[]>& epsilonCol, int n){
    std::cout << " second phase begin " << std::endl;
    int d = -1; 
    int r; 
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(epsilonCol[j] == 0 && epsilonRow[i] == 0){
                r = u[i] + v[j] - R[i*n + j] ; 
                if(d == -1 || d > r ) {
                    d = r; 
                } 
            }
        } 
    }
    if(d == -1 || d == 0){
        printSquaredMatrix("The solution is : ", Q, n*n);  
        return ;
    } else {
        int m = -1;
        int mz = -1 ; 
        bool firstRow = true;
        for(int i =0 ; i < n; i++) {
            if(epsilonRow[i] == 0){
                if(m == -1 || (u[i] < m  && u[i] > 0)){
                    m = u[i];  
                }  
                if(u[i] == 0 && firstRow){
                    firstRow = false; 
                    for(int j = 0 ; j < n ; j++){
                        if(epsilonCol[j] == 0) {
                            if(mz == -1 || (v[j] < mz)) {
                                mz = v[j]; 
                            }
                        }
                    }
                    if(d < mz) {
                        mz = d;  
                    }
                }
            }
        }
        if(m > d) 
            m = d; 
        for(int i=0 ; i < n ; i++){
            if(epsilonRow[i] == 0){
                u[i] -= m;
            } else {
                u[i] += mz;
            }
            if(epsilonCol[i] == 1){
                v[i] += m; 
            } else {
                v[i] -= mz;
            }

        }
        printSquaredMatrix("u after second phase ", u, n);
        printSquaredMatrix("v after second phase ", v, n);
        firstPhase(R, Q, u, v, ik, jk, epsilonRow, n, OTHERSTAGE);

    }
}
/**
 * Q: pointer to the associated Q matrix of the paper.  
 * epsilon: pointer to the record of essential rows.  
 * n: size of rows/cols. (nxn) matrix.
 */
void getEssentialCol(const std::shared_ptr<int[]>& Q, const std::shared_ptr<int[]> &epsilon, const std::shared_ptr<int[]> &result, int n){
    for(int i =0 ; i < n ; i++){
        if(epsilon[i] == 0){
            for(int j = 0; j < n ; j++) {
                if(Q[i*n+j] == 2) {
                    result[j] = 1; 
                }
            }
        }
    }
}

int main(){
    int n; 
    std::cin >> n; 
    n = 4;
    const std::shared_ptr<int[]> flaged(new int[n*n]);
    const std::shared_ptr<int[]> R (new int[n*n]);
    const std::shared_ptr<int[]> Q (new int[n*n]);
    const std::shared_ptr<int[]> u (new int[n]);
    const std::shared_ptr<int[]> v (new int[n]);
    const std::shared_ptr<int[]> epsilonRow (new int[n]);
    fillCustom(epsilonRow, n, 0);
    const std::shared_ptr<int[]> epsilonCol (new int[n]);
    fillCustom(epsilonCol, n, 0);


    fillCustom(R, n*n,  n);
    R[0]=8; 
    R[1]=7; 
    R[2]=9; 
    R[3]=9; 
    R[4]=5; 
    R[5]=2; 
    R[6]=7;
    R[7]=8;
    R[8]=6;
    R[9]=1;
    R[10]=4;
    R[11]=9;
    R[12]=2;
    R[13]=3;
    R[14]=2;
    R[15]=6;
    fillCustom(flaged, n*n, 0);

    printSquaredMatrix("R", R, n*n);
    //printSquaredMatrix(flaged, 3*3);

    const std::shared_ptr<int[]> ik(new int[n]);
    const std::shared_ptr<int[]> jk(new int[n]);
    fillCustom(ik, n, -1);
    fillCustom(jk, n, -1);

    firstPhase(R, Q, u, v, ik, jk,  epsilonRow, n, FIRSTSTAGE);
    printSquaredMatrix("epsilonRow", epsilonRow, n);

    getEssentialCol(Q, epsilonRow, epsilonCol, n);
    printSquaredMatrix("epsilonCol", epsilonCol, n);

    //void secondPhase(std::unique<int[]>& u, std::unique<int[]>& v, const std::shared_ptr<int[]>& epsilonRow, const std::shared_ptr<int[]>& epsilonCol, int n){
    printSquaredMatrix("u", u, n);
    printSquaredMatrix("v", v, n);

    secondPhase(R, Q, u, v, ik, jk,  epsilonRow, epsilonCol, n);
    return 0;
}
