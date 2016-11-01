#include <stdio.h>
#include <math.h>


float polynom ( float A1 , float A2 , float A3 , float A4 , float A5 , float A6 , float A7 , float z , float n) {
    float F;
    F = A1 * powf ( z , 2*n ) + A2 * powf ( z , n+2 ) + A3 * powf ( z , n+1 ) + A4 * powf ( z , n ) + A5 * powf (z , 2 ) + A6 * z + A7; // искомое уравнение
    return F;
}

void solution ( float Q1 , float S1 , float A1 , float A2 , float A3 , float A4 , float A5 , float A6 , float A7 , float n , float delta , float U3 , float U0 , float r0 , float P0 , float P3) {
    int i = 0;
    float K; float z1 ; float z2; float z3; float U; float D0;
    
    
    K = fabsf((S1 - Q1))/2 + Q1;
    
    
    while ( i == 0 ) {
        z1 = polynom ( A1 , A2 , A3 , A4 , A5 , A6, A7 , Q1 , n );
        z2 = polynom ( A1 , A2 , A3 , A4 , A5 , A6, A7 , K , n );
        z3 = polynom ( A1 , A2 , A3 , A4 , A5 , A6, A7 , S1 , n );
        
        
        if ( z1 > 0 ) {
            if ( z2 < 0 ) {
                if ( fabsf ( z2 ) <= 0.0001 ) {
                    i = 1;
                    printf ( " solution = %.20f\n" ,K );
                    U = U3 + delta*(1-K );
                    D0 = (r0*powf(U0 , 2) + P0 -powf(K , n )*P3 - r0*U0*U)/(r0*U0-r0*U);
                    printf ( " D0 = %.10f cm/c\n" , D0);
                    
                } else {
                    
                    S1 = K;
                    K = Q1 + fabsf( K - Q1)/2;
                    
                }
            }
        }
        if ( z1 < 0 ) {
            if ( z2 > 0 ) {
                if ( fabsf ( z2 ) <= 0.0001 ) {
                    i = 1;
                    printf ( " solution = %.20f\n" , K  );
                    U = U3 + delta*(1-K);
                    D0 = (r0*powf(U0 , 2) + P0 -powf(K , n )*P3 - r0*U0*U)/(r0*U0-r0*U);
                    printf ( " D0 = %.10f cm/c\n" , D0);
                } else {
                    
                    S1 = K;
                    K = Q1 + fabsf( K - Q1)/2;
                    
                }
            }
        }
        if ( z2 > 0 ) {
            if ( z3 < 0 ) {
                if ( fabsf ( z2 ) <= 0.0001 ) {
                    i = 1;
                    printf ( " solution = %.20f\n" , K  );
                    U = U3 + delta*(1-K);
                    D0 = (r0*powf(U0 , 2) + P0 -powf(K , n )*P3 - r0*U0*U)/(r0*U0-r0*U);
                    printf ( " D0 = %.10f cm/c\n" , D0);
                } else {
                    
                    Q1 = K;
                    K = K + fabsf( K - S1)/2;
                    
                }
            }
        }
        if ( z2 < 0 ) {
            if ( z3 > 0 ) {
                if ( fabsf ( z2 ) <= 0.0001 ) {
                    i = 1;
                    printf ( " solution = %.20f\n" , K  );
                    U = U3 + delta*(1-K);
                    D0 = (r0*powf(U0 , 2) + P0 -powf(K , n )*P3 - r0*U0*U)/(r0*U0-r0*U);
                    printf ( " D0 = %.10f cm/c\n" , D0);
                } else {
                    
                    Q1 = K;
                    K = K + fabsf( K - S1)/2;
                }
            }
        }
    }
}

void local ( float A1 , float A2 , float A3 , float A4 , float A5 , float A6 , float A7 , float n , float delta, float U3 , float U0 , float r0 , float P0 , float P3 ) {
    
    float A; float B; int i; float M[7]; int N = 10000000; float L ; float R; float S; float Q; float z1 ; float z2; float Q1; float S1;		//  локализируем корни
    
    M[0] = A1; M[1] = A2; M[2] = A3; M[3] = A4; M[4] = A5; M[5] = A6; M[6] = A7;
    
    A = 0; B = 0;
    
    for ( i = 1 ; i < 6 ; i++ ) {     // определяем максимальную область
        if (A < fabsf ( M[i] )) {
            A = fabsf(M[i]);
        }
    }
    
    for (i=0 ; i<5 ; i++ ) {
        if ( B < fabsf ( M[i] )) {
            B = fabsf(M[i]);
        }
    }
    
    printf("\n");
    printf ( "max A = %f , max B = %f\n " , A , B);
    
    L =  (fabsf(M[6]) / (B + fabsf(M[6]) )) ;
    R = ( 1 + A / fabs ( M[0] ) ) ;
    
    printf("\n");
    printf( " L = %f , R = %f\n" , L , R );
    
    S = ((R-L)/N );			// задаем шаг
    i = 1; Q = L;
    printf ( "S = %f\n" , S );
    
    printf ("\n");
    
    printf ( " область локализации корней: \n " );
    while ( Q+S < R) {
        z1 = polynom ( A1 , A2 , A3 , A4 , A5 , A6 , A7 , Q , n );
        z2 = polynom ( A1 ,  A2 ,  A3 , A4 ,  A5 , A6 , A7 , Q + S , n );	// поиск области при положительных значениях
        
        if ( z1 > 0 ) {
            if ( z2 < 0 ) {
                printf ( " %d = [ %.20f , %.20f ]\n " ,i, Q , Q + S );
                printf ( "корень:  %f ", Q + S/2);
                Q1 = Q;
                S1 = Q + S;
                //solution(Q1 , S1 , A1 , A2 , A3 , A4 , A5 , A6 , A7 , n , delta , U3 ,  U0, r0, P0, P3);
                
                i++;
            }
        }
        if  (z1 < 0 ) {
            if  ( z2 > 0 )  {
                printf ( " %d = [ %.20f , %.20f ]\n " ,i, Q , Q + S );
                printf ( "корень: %f ", Q + S/2);
                Q1 = Q;
                S1 = Q + S;
                //solution(Q1 , S1 , A1 , A2 , A3 , A4 , A5 , A6 , A7 , n , delta , U3 ,  U0, r0, P0, P3);
                
                i++;
            }
        }
        Q = Q + S;
    }
    
    
}

int main() {
    float j0 , r0 , U0 , P0 , j3 , C3 , U3 , P3; // conditions at the beginning
    float a0 , n ,m , v , x ; // unknown variables + z and P1
    float A1 , A2 , A3 , A4 , A5 , A6 , A7; // parts of polynom
    //printf ( "vvedite j0 , r0 , U0 , P0 , j3 , C3 , U3 , P3\n" );
    float delta;
    
    j0 = 1.4 ; r0 = 0.0001694 ; U0 = 0.001 ; P0 = 1013000 ; j3 = 1.4 ; C3 = 36537 ; U3 = 12290 ; P3 = 1676800;
    
    
    //scanf ( "%f , %f , %f , %f , %f , %f , %f , %f"  , &j0 , &r0 , &U0 , &P0 , &j3 , &C3 , &U3 , &P3 ); // For other variants of conditions
    
    
    
    a0 = ( j0 + 1 )/( j0 - 1 );
    n = 2*j3/( j3 - 1 );
    m = ( U3 - U0 )*   sqrtf(( j0 - 1 ) * r0 / (2 * P0 ))   ;    // NEED TO CHEK PHYSICAL DIMENTIONS OF VARIABLES!!!!!!
    v = (2*C3/(j3 - 1))* sqrtf ( r0*(j0-1)/(2*P0));
    x = P3/P0;
    delta = (2*C3/(j3 - 1));
    printf ( "a0 = %f ,n = %f ,m = %f ,v = %f ,x = %f\n" , a0 , n , m , v , x );
    printf( "\n" );
    A1 = powf ( x , 2 );
    A2 = (-1)*a0*powf( v , 2 )*x;
    A3 = 2*a0*v*(m+v)*x;
    A4 = (-1)*(2+ powf ( m+v , 2 )*a0)*x;
    A5 = (-1)*powf( v , 2 );
    A6 = 2*v*(m+v);
    A7 = (-1)*powf(m+v, 2) +1;
    printf ( "A1 = %f ,A2 = %f ,A3 = %f ,A4 = %f ,A5 = %f ,A6 = %f ,A7 = %f\n" , A1 , A2 , A3 , A4 , A5 , A6 , A7 );
    local ( A1 , A2 , A3, A4 , A5, A6, A7 , n , delta , U3 , U0, r0, P0, P3);
    
    return 0;
}

