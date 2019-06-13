#include <math.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
// structure permettant de definir Fp2
typedef struct{
    mpz_t x;
    mpz_t y;
   }comp;
typedef comp comp_t[1];
// Entetes des fonctions
void comp_init(comp_t);
void add_ext(comp_t ,comp_t ,comp_t ,mpz_t );
void mul_ext(comp_t,comp_t,comp_t,mpz_t);
void sub_ext(comp_t,comp_t,comp_t,mpz_t);
void inv_ext(comp_t,comp_t,mpz_t);
void set_ext(comp_t,comp_t);
void opposer(comp_t,comp_t,mpz_t);
void racine(mpz_t,mpz_t,mpz_t);
void root_ext(comp_t,comp_t,mpz_t);
void add_ext_ui(comp_t,comp_t,int,mpz_t);
void set_ext_ui(comp_t,int);
void add_ext_ui(comp_t,comp_t,int,mpz_t);
void sub_ext_ui(comp_t,comp_t,int,mpz_t);
void xDBL(comp_t,comp_t,comp_t,comp_t,comp_t,comp_t ,mpz_t);
void double_and_add(comp_t ,comp_t ,char*,comp_t,comp_t,comp_t,comp_t,mpz_t);
void xDBLe(comp_t , comp_t ,comp_t ,comp_t ,comp_t , comp_t , int ,mpz_t );
void iso_2_e(comp_t,comp_t,comp_t,comp_t,comp_t,comp_t ,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,
               comp_t,comp_t,int,mpz_t);
void xTPL(comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,mpz_t);
void xTPLe(comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,int,mpz_t);
void xADD(comp_t,comp_t, comp_t,comp_t,comp_t,comp_t,comp_t,comp_t ,mpz_t);
void get_xR(comp_t ,comp_t,comp_t ,comp_t,comp_t,comp_t ,comp_t,mpz_t);
void j_inv(comp_t,comp_t,comp_t,mpz_t);
void curve_4_iso(comp_t,comp_t,comp_t,comp_t,comp_t,mpz_t);
void point_4_iso(comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,mpz_t);
void get_A(comp_t,comp_t,comp_t,comp_t,mpz_t);
void get_yP_yQ_a_b(comp_t ,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t, mpz_t);
void curve_3_iso(comp_t,comp_t,comp_t,comp_t,comp_t,mpz_t);
void point_3_iso(comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,mpz_t);
void iso_3_e(comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,
             comp_t,comp_t,int,mpz_t);

void isogen2(comp_t ,comp_t,comp_t,mpz_t,int,int,comp_t,comp_t,comp_t,comp_t,
                         comp_t,comp_t,comp_t,comp_t ,comp_t,comp_t,mpz_t);

void isogen3(comp_t,comp_t,comp_t, mpz_t,int,int,comp_t ,comp_t ,comp_t,comp_t,
                comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,mpz_t);

void isoex2(comp_t,mpz_t,comp_t,comp_t,comp_t,int,int,comp_t,comp_t,mpz_t);

void isoex3(comp_t,mpz_t,comp_t,comp_t,comp_t,int,int,comp_t,comp_t,mpz_t);

int main(){


//déclaration des variables
comp_t a,b;
mpz_init(a->x);mpz_init(a->y);mpz_init(b->x);mpz_init(b->y);
   //les coéfficients de Montgoméry de la courbe
int e2,e3; //les exposants publiques de 2 et 3 respectivement
mpz_t p;
mpz_init(p);//le nombre premier
mpz_t x,y;
mpz_init(x);mpz_init(y);

mpz_t sk2,sk3;//les entiers aléatoires secrets d'Alice et Bob respectivement
mpz_init(sk2);mpz_init(sk3);
gmp_randstate_t state;
gmp_randinit_default (state);
comp_t xP2,yP2,xQ2,yQ2;//les coordonnés des points pubiques d'Alice
comp_t xP3,yP3,xQ3,yQ3;//les coordonnées des points publiques de Bob
mpz_init(xP2->x);mpz_init(xP2->y);mpz_init(yP2->x);mpz_init(yP2->y);
mpz_init(xQ2->x);mpz_init(xQ2->y);mpz_init(yQ2->x);mpz_init(yQ2->y);
mpz_init(xP3->x);mpz_init(xP3->y);mpz_init(yP3->x);mpz_init(yP3->y);
mpz_init(xQ3->x);mpz_init(xQ3->y);mpz_init(yQ3->x);mpz_init(yQ3->y);
comp_t xP2_prime,xQ2_prime,xR2_prime;//Public key pk3
comp_t xP3_prime,xQ3_prime,xR3_prime;//Public key pk2
comp_init(xP2_prime);comp_init(xQ2_prime);comp_init(xR2_prime);
comp_init(xP3_prime);comp_init(xQ3_prime);comp_init(xR3_prime);
printf("\n******************************************************************************");
printf("           :Entrée des paramètres publiques  :\n");
printf("\n******************************************************************************");

printf("\nEntrer l'exposant e2 de 2: ");
scanf("%d",&e2);
printf("\nEntrer l'exposant e3 de 3 : ");
scanf("%d",&e3);

printf("\nEntrer les coordonnées du coefficient a de la courbe:");
gmp_scanf("%Zd %Zd",&x,&y);
mpz_set(a->x,x);mpz_set(a->y,y);

printf("\nEntrer les coordonnées du coefficient b de la courbe:");
gmp_scanf("%Zd %Zd",&x,&y);
mpz_set(b->x,x);mpz_set(b->y,y);

printf("\nEntrer l'abscisse du point P2 choisit par Alice:");
gmp_scanf("%Zd %Zd",&x,&y);
mpz_set(xP2->x,x);mpz_set(xP2->y,y);
printf("\nEntrer l'ordonnée du point P2 choisit par Alice:");
gmp_scanf("%Zd %Zd",&x,&y);
mpz_set(yP2->x,x);mpz_set(yP2->y,y);
printf("\nEntrer l'abscisse du point Q2 choisit par Alice:");
gmp_scanf("%Zd %Zd",&x,&y);
mpz_set(xQ2->x,x);mpz_set(xQ2->y,y);
printf("\nEntrer l'ordonnée du point Q2 choisit par Alice:");
gmp_scanf("%Zd %Zd",&x,&y);
mpz_set(yQ2->x,x);mpz_set(yQ2->y,y);

printf("\nEntrer l'abscisse  du point P3 choisit par Bob :");
gmp_scanf("%Zd %Zd",&x,&y);
mpz_set(xP3->x,x);mpz_set(xP3->y,y);
printf("\nEntrer l'ordonnée  du point P3 choisit par Bob :");
gmp_scanf("%Zd %Zd",&x,&y);
mpz_set(yP3->x,x);mpz_set(yP3->y,y);
printf("\nEntrer l'abscisse  du point Q3 choisit par Bob :");
gmp_scanf("%Zd %Zd",&x,&y);
mpz_set(xQ3->x,x);mpz_set(xQ3->y,y);
printf("\nEntrer l'ordonnée du point Q3 choisit par Bob:");
gmp_scanf("%Zd %Zd",&x,&y);
mpz_set(yQ3->x,x);mpz_set(yQ3->y,y);

mpz_clear(x);mpz_clear(y);
mpz_t puiss2,puiss3,prod,var;
mpz_init(puiss2);mpz_init(puiss3);mpz_init(prod);mpz_init(var);
mpz_set_ui(var,2);
mpz_pow_ui(puiss2,var,e2);
mpz_set_ui(var,3);
mpz_pow_ui(puiss3,var,e3);
mpz_mul(prod,puiss2,puiss3);
mpz_mul_ui(prod,prod,11);
mpz_set_ui(var,1);
mpz_sub(p,prod,var);

printf("\n******************************************************************************\n");
printf("\n Générations de la clé publique d'Alice :\n");
printf("\n******************************************************************************\n");
mpz_urandomm(sk2,state,puiss2);

/*mul_ext(yP2,xP2,xP2,p);
mul_ext(yP2,yP2,xP2,p);
add_ext(yP2,yP2,xP2,p);
racine(yP2->x,yP2->x,p);
mpz_set_ui(yP2->y,0);
xTPLe(xP2,yP2,xP2,yP2,a,b,e3,p);
set_ext(xQ2,xP2); opposer(xQ2,xP2,p);
mpz_set(yQ2->y,yP2->x);mpz_add(yQ2->x,yP2->y,yP2->y);mpz_sub(yQ2->x,yQ2->y,yQ2->x);*/

isogen2(xP3_prime,xQ3_prime,xR3_prime,sk2,e2,e3,xP2,yP2,xQ2,yQ2,xP3,yP3,xQ3,yQ3,a,b,p);

gmp_printf ("la clé publique d'Alice est pk2=(%Zd + %Zdi/%Zd + %Zdi/%Zd + %Zdi)\n:",
            xP3_prime->x,xP3_prime->y,xQ3_prime->x,xQ3_prime->y,xR3_prime->x,xR3_prime->y);



printf("\n******************************************************************************\n");
printf(" \nGénérations de la clé publique de Bob:\n");
printf("\n******************************************************************************\n");

mpz_urandomm(sk3,state,puiss3);
/*mul_ext(yP3,xP3,xP3,p);
mul_ext(yP3,yP3,xP3,p);
add_ext(yP3,yP3,xP3,p);
racine(yP3->x,yP3->x,p);
mpz_set_ui(yP3->y,0);

xDBLe(xP3,yP3,xP3,yP3,a,b,e2,p);
set_ext(xQ3,xP3); opposer(xQ3,xQ3,p);
mpz_set(yQ3->y,yP3->x);mpz_add(yQ3->x,yP3->y,yP3->y);mpz_sub(yQ3->x,yP3->y,yQ3->x);*/


isogen3(xP2_prime,xQ2_prime,xR2_prime,sk3,e2,e3,xP2,yP2,xQ2,yQ2,xP3,yP3,xQ3,yQ3,a,b,p);


gmp_printf ("la clé publique de Bob  est pk3=(%Zd + %Zdi/%Zd + %Zdi/%Zd + %Zdi \n)",
            xP2_prime->x,xP2_prime->y,xQ2_prime->x,xQ2_prime->y,xR2_prime->x,xR2_prime->y);

printf("\n******************************************************************************\n");
printf(" \n Calcul du secret partagé par Alice   :\n");
printf("\n******************************************************************************\n");

comp_t j2;//le j_invariant de la courbe secrète d'Alice
comp_init(j2);
isoex2(j2,sk2,xP3_prime,xQ3_prime,xR3_prime,e2,e3,a,b,p);

gmp_printf ("la clé privée d'Alice est j2=%Zd + %Zdi \n",j2->x,j2->y);

printf("\n******************************************************************************\n");
printf("\nCalcul du secret partagé  par Bob   :\n");
printf("\n******************************************************************************\n");

comp_t j3;//le j_invariant de la courbe secrète d'Alice
comp_init(j3);
isoex3(j3,sk3,xP2_prime,xQ2_prime,xR2_prime,e2,e3,a,b,p);

gmp_printf ("la clé privée d'Alice est j3=%Zd + i%Zd \n",j3->x,j3->y);

mpz_mod(j2->x,j2->x,p);mpz_mod(j2->y,j2->y,p);
mpz_mod(j3->x,j3->x,p);mpz_mod(j3->y,j3->y,p);
if(mpz_cmp(j2->x,j3->x)==0 && mpz_cmp(j2->y,j3->y)==0)
    printf("Echange réussie\n");
else
    printf("echec:\n");
return 0;
}
//initialisation d
void comp_init(comp_t C){
   mpz_init(C->x);
   mpz_init(C->y);
}

//addition dans Fp2
void add_ext(comp_t res, comp_t a,comp_t b,mpz_t p)
 {
   mpz_add(res->x,a->x,b->y);
   mpz_add(res->y,a->x,b->y);
   mpz_mod(res->x,res->x,p);
   mpz_mod(res->y,res->y,p);
}
//soustraction dans fp2
void sub_ext(comp_t res,comp_t a,comp_t b,mpz_t p)
  {
    mpz_sub(res->x,a->x,b->x);
    mpz_sub(res->y,a->y,b->y);
    mpz_mod(res->x,res->x,p);
    mpz_mod(res->y,res->y,p);
  }
//multiplication dans Fp2
  void mul_ext(comp_t res,comp_t t1 ,comp_t t2,mpz_t p)
   { mpz_t inter1,inter2;
     mpz_init(inter1);mpz_init(inter2);
     mpz_mul(res->x,t1->x,t2->x);
     mpz_mul(inter1,t1->y,t2->y);
     mpz_sub(res->x,res->x,inter1);
     mpz_mul(res->y,t1->x,t2->y);
     mpz_mul(inter2,t2->x,t1->y);
     mpz_add(res->y,res->y,inter2);
     mpz_mod(res->x,res->x,p);
     mpz_mod(res->y,res->y,p);
   }
// l'opposer d'un element de Fp2: pour a €fp2,opposer calcul -a

   void opposer(comp_t res,comp_t a,mpz_t p)
     {
       add_ext(res,a,a,p);
       sub_ext(res,a,res,p);
     }
     void set_ext(comp_t res,comp_t a)
       {
         mpz_set(res->x,a->x);
         mpz_set(res->y,a->y);
         return;
       }

       void set_ext_ui(comp_t res,int a)
         {  comp_t b;
           comp_init(b);
           mpz_set_ui(b->x,a);
           mpz_set_ui(b->y,0);
           set_ext(res,b);
        }

// addition d'un élément de Fp2 et d'un element de type mpz
      void add_ext_ui(comp_t res,comp_t a,int b,mpz_t p)
         {
           comp_t c;
           comp_init(c);
           mpz_set_ui(c->x,b);
           mpz_set_ui(c->y,0);
           add_ext(res,a,c,p);
         }
//soustraction d'un élément de Fp2 et d'un element de type mpz
     void sub_ext_ui(comp_t res,comp_t a,int b,mpz_t p)
        {
          comp_t c;
          comp_init(c);
          mpz_set_ui(c->x,b);
          mpz_set_ui(c->y,0);
          sub_ext(res,a,c,p);
        }


// inversion d'un élément de Fp2
        void inv_ext(comp_t res,comp_t t,mpz_t p)
         { mpz_t a,b,module;
           mpz_init(a);mpz_init(b);mpz_init(module);
           mpz_mul(a,t->x,t->x);
           mpz_mul(b,t->y,t->y);
           mpz_add(module,a,b);
           mpz_invert(module,module,p);
           mpz_mul(res->x,module,t->x);
           mpz_mul(res->y,module,t->y);
           mpz_mod(res->x,res->x,p);
           mpz_mod(res->y,res->y,p);
       }
//calcul de la racine carré d'un élément de type mpz
                  void racine(mpz_t res,mpz_t a,mpz_t p){
                       mpz_t p1,exp;
                       mpz_init(p1);mpz_init(exp);
                       mpz_add_ui(p1,p,1);
                       mpz_divexact_ui(exp,p1,4);
                       mpz_powm(res,a,exp,p);
                  }
//calcul de la racine carré d'un élément de Fp2
       void root_ext(comp_t res,comp_t t1,mpz_t p)
         {
           mpz_t a,b,prod,b1,t2,deux,delta;
           mpz_init(a);mpz_init(b);mpz_init(prod);mpz_init(delta);
           mpz_init(b1);mpz_init(t2);mpz_init(deux);
           mpz_set(a,t1->x);
           mpz_set(b,t1->y);
           mpz_set_ui(deux,2);
           mpz_mul(delta,a,a);
           mpz_mul(b1,b,b);
           mpz_sub(delta,delta,b1);
           racine(delta,delta,p);
           mpz_add(delta,a,delta);
           mpz_invert(deux,deux,p);
           mpz_mul(res->x,delta,deux);
           racine(res->x,res->x,p);

             mpz_sub(t2,a,delta);
             mpz_mul_ui(res->y,res->x,2);
             mpz_invert(res->y,res->y,p);
             mpz_mul(res->y,res->y,a);
          }

// calcul de 2P pour P un point de la courbe elliptique
       void xDBL(comp_t x2P,comp_t y2P,comp_t xP,comp_t yP,comp_t a,comp_t b,mpz_t p)
       {
         comp_t t0,t1,t2;
         comp_init(t0);comp_init(t1);comp_init(t2);
         mul_ext(t0,xP,xP,p);
         add_ext(t1,t0,t0,p);
         add_ext(t0,t0,t1,p);
         mul_ext(t1,a,xP,p);
         add_ext(t1,t1,t1,p);
         add_ext(t0,t0,t1,p);
         add_ext(t0,t0,t2,p);
         mul_ext(t1,b,yP,p);
         add_ext(t1,t1,t1,p);
         inv_ext(t1,t1,p);
         mul_ext(t0,t0,t1,p);
         mul_ext(t1,t0,t0,p);
         mul_ext(t2,b,t1,p);
         sub_ext(t2,t2,a,p);
         sub_ext(t2,t2,xP,p);
         sub_ext(t2,t2,xP,p);
         mul_ext(t1,t0,t1,p);
         mul_ext(t1,b,t1,p);
         add_ext(t1,t1,yP,p);
         add_ext(y2P,xP,xP,p);
         add_ext(y2P,y2P,xP,p);
         add_ext(y2P,y2P,a,p);
         mul_ext(y2P,y2P,t0,p);
         sub_ext(y2P,y2P,t1,p);
         set_ext(x2P,t2);
       }
//calcul de kP :ici le tableau m represente les bits de k
         void double_and_add(comp_t xmP,comp_t ymP,char* m,comp_t x,comp_t y,comp_t a,comp_t b,
                           mpz_t p)
           { int l=strlen(m);
             int i;
             comp_t x0,y0;
             comp_init(x0);comp_init(y0);
             set_ext_ui(x0,0);
             set_ext_ui(y0,0);
             for (i=l-1;i<=0;i++){
                      xDBL(x0,y0,x0,y0,a,b,p);
                      if(m[i]==1)
                         xADD(x0,y0,x0,y0,x,y,a,b,p);

                   }
               }
//calcul de (2^e)P,P etant un point de la courbe elliptique
         void xDBLe(comp_t xR, comp_t yR,comp_t xP,comp_t yP,comp_t a, comp_t b, int e,mpz_t p)
             {
               comp_t x ,y,res_x,res_y;
               comp_init(x);comp_init(y);comp_init(res_x);comp_init(res_y);
               set_ext(x,xP);
               set_ext(y,yP);
               int i=0;
               while(i<e){
                          xDBL(res_x,res_y,x,y,a,b,p);
                          set_ext(x,res_x);
                          set_ext(y,res_y);
                          i++;

                         }
               set_ext(xR,x);
               set_ext(yR,y);
               return;
           }
//calcul l'image de la courbe de départ par une isogénie de degré 2
       void iso_2_e(comp_t a_res,comp_t b_res,comp_t res_xP1,comp_t res_yP1,comp_t res_xP2,
                   comp_t res_yP2,comp_t xS,comp_t yS,comp_t xP1,comp_t yP1,comp_t xP2,comp_t yP2,
                    comp_t a,comp_t b,int e2,mpz_t p)
           {
             int e=e2-2;
             int i=e;
             comp_t xT,yT;
             comp_init(xT);comp_init(yT);
             while(i>0){
                    xDBLe(xT,yT,xS,yS,a,b,e,p);
                    curve_4_iso(a_res,b_res,xT,a,b,p);
                    point_4_iso(xS,yS,xS,yS,xT,a,b,p);
                    i=i-2;
                      }

                   point_4_iso(res_xP1,res_yP1,xP1,yP1,xT,a,b,p);
                   point_4_iso(res_xP2,res_yP2,xP2,yP2,xT,a,b,p);
             return;
           }
//calcul de 3P ,P etant un point de la courbe elliptique

             void xTPL(comp_t x3P,comp_t y3P,comp_t xP,comp_t yP,comp_t a,comp_t b,mpz_t p)
                   {
                     comp_t x2P,y2P;
                     comp_init(x2P);comp_init(y2P);
                     xDBL(x2P,y2P,xP,yP,a,b,p);
                     xADD(x3P,y3P,xP,yP,x2P,x2P,a,b,p);
                   }
//calcul de (3^e)P , P etant un point de la courbe
             void xTPLe(comp_t x3P,comp_t y3P,comp_t xP,comp_t yP,comp_t a,comp_t b,int e,mpz_t p)
                {
                 comp_t x,y,res_x,res_y;
                 comp_init(x);comp_init(y);comp_init(res_x);comp_init(res_y);
                 set_ext(x,xP);
                 set_ext(y,yP);
                 int i=0;
                 while(i<e){
                            xTPL(res_x,res_y,x,y,a,b,p);
                            set_ext(x,res_x);
                            set_ext(y,res_y);
                            i++;
                           }
                 set_ext(x3P,x);
                 set_ext(y3P,y);
               return;
             }
// calcul de 2p
             void xADD(comp_t xPQ,comp_t yPQ, comp_t xP,comp_t yP,comp_t xQ,comp_t yQ,comp_t a,comp_t b,mpz_t p)
             {
               if((mpz_cmp(xP->x,xQ->x)==0)&&(mpz_cmp(xP->y,xQ->y)==0)&&(mpz_cmp(yP->x,yQ->x)==0)&&
               (mpz_cmp(yP->y,yQ->y)==0)){
                                       xDBL(xPQ,yPQ,xP,yP,a,b,p);}
               else { comp_t t0,t1,t2;
                      comp_init(t0);comp_init(t1);comp_init(t2);
                      mul_ext(t0,b,t0,p);
                      add_ext(t0,t0,yP,p);
                      sub_ext(t0,t2,t0,p);
                      mul_ext(t1,b,t1,p);
                      sub_ext(t1,t1,a,p);
                      sub_ext(t1,t1,xP,p);
                      sub_ext(xPQ,t1,xQ,p);
                      set_ext(yPQ,t0);
                    }
               return;
             }
    //calcul l'abscisse  du point R=  P-Q
             void get_xR(comp_t xR,comp_t a,comp_t b,comp_t xP,comp_t yP,comp_t xQ,comp_t yQ,mpz_t p)
                 { comp_t yR;
                   comp_init(yR);
                   opposer(yQ,yQ,p);
                   xADD(xR,yR,xP,yP,xQ,yQ,a,b,p);
                   return;
                 }
  //calcul du j-invariant d'une courbe elliptique
                 void j_inv(comp_t j,comp_t a,comp_t b,mpz_t p)
                       {
                         comp_t t0,t1;
                         comp_init(t0);comp_init(t1);
                         mul_ext(t0,a,a,p);
                         set_ext_ui(j,3);
                         sub_ext(j,t0,j,p);
                         mul_ext(t1,j,j,p);
                         mul_ext(j,j,t1,p);
                         add_ext(j,j,j,p);
                         add_ext(j,j,j,p);
                         add_ext(j,j,j,p);
                         add_ext(j,j,j,p);
                         add_ext(j,j,j,p);
                         add_ext(j,j,j,p);
                         add_ext(j,j,j,p);
                         add_ext(j,j,j,p);
                         set_ext_ui(t1,4);
                         sub_ext(t0,t0,t1,p);
                         inv_ext(t0,t0,p);
                         mul_ext(j,j,t0,p);
                         return;

                       }
                       //calcul de l'image de la courbe elliptique par l'isogenie de de degré 4
                       void curve_4_iso(comp_t a_image,comp_t b_image,comp_t xP4,comp_t a,comp_t b,mpz_t p)
                          {
                           comp_t t1,t2;
                           comp_init(t1);comp_init(t2);
                           mul_ext(t1,xP4,xP4,p);
                           mul_ext(a_image,t1,t1,p);
                           add_ext(a_image,a_image,a_image,p);
                           add_ext(a_image,a_image,a_image,p);

                           set_ext_ui(t2,2);
                           sub_ext(a_image,a_image,t2,p);
                           mul_ext(t1,xP4,t1,p);
                           add_ext(t1,t1,xP4,p);
                           mul_ext(t1,t1,b,p);
                           inv_ext(t2,t2,p);
                           opposer(t2,t2,p);
                           mul_ext(b_image,t2,t1,p);
                           return;
                         }
//calcul de l'image d'un point par l'isogénie de degré 4
                       void point_4_iso(comp_t xQres,comp_t yQres,comp_t xQ,comp_t yQ,comp_t xP4,
                                          comp_t a, comp_t b,mpz_t p)
                         {
                           comp_t t1,t2,t3,t4,t5;
                           comp_init(t1);comp_init(t2);comp_init(t3);comp_init(t4);comp_init(t5);
                           mul_ext(t1,xQ,xQ,p);
                           mul_ext(t2,t1,t1,p);
                           mul_ext(t3,xP4,xP4,p);
                           mul_ext(t4,t2,t3,p);
                           add_ext(t2,t2,t4,p);
                           mul_ext(t4,t1,t3,p);
                           add_ext(t4,t4,t4,p);
                           add_ext(t5,t4,t4,p);
                           add_ext(t5,t5,t5,p);
                           add_ext(t4,t4,t5,p);
                           add_ext(t2,t2,t4,p);
                           mul_ext(t4,t3,t3,p);
                           mul_ext(t5,t1,t4,p);
                           add_ext(t5,t5,t5,p);
                           add_ext(t2,t2,t5,p);
                           mul_ext(t1,t1,xQ,p);
                           mul_ext(t4,xP4,t3,p);
                           mul_ext(t5,t1,t4,p);
                           add_ext(t5,t5,t5,p);
                           add_ext(t5,t5,t5,p);
                           sub_ext(t2,t2,t5,p);
                           mul_ext(t1,t1,xP4,p);
                           add_ext(t1,t1,t1,p);
                           add_ext(t1,t1,t1,p);
                           sub_ext(t1,t2,t1,p);
                           mul_ext(t2,xQ,t4,p);
                           add_ext(t2,t2,t2,p);
                           add_ext(t2,t2,t2,p);
                           sub_ext(t1,t1,t2,p);
                           add_ext(t1,t1,t3,p);
                           add_ext_ui(t1,t1,1,p);
                           mul_ext(t2,xQ,xP4,p);
                           sub_ext_ui(t4,t2,1,p);
                           add_ext(t2,t2,t2,p);
                           add_ext(t5,t2,t2,p);
                           sub_ext(t1,t1,t5,p);
                           mul_ext(t1,t1,t4,p);
                           mul_ext(t1,t1,t3,p);
                           mul_ext(t1,yQ,t1,p);
                           add_ext(t1,t1,t1,p);
                           opposer(yQres,t1,p);
                           sub_ext(t2,t2,t3,p);
                           sub_ext_ui(t1,t2,1,p);
                           sub_ext(t2,xQ,xP4,p);
                           mul_ext(t1,t1,t2,p);
                           mul_ext(t5,t1,t1,p);
                           mul_ext(t5,t2,t5,p);
                           inv_ext(t5,t5,p);
                           mul_ext(yQres,yQres,t5,p);
                           mul_ext(t1,t1,t2,p);
                           inv_ext(t1,t1,p);
                           mul_ext(t4,t4,t4,p);
                           mul_ext(t1,t1,t4,p);
                           mul_ext(t1,xQ,t1,p);
                           mul_ext(t2,xQ,t3,p);
                           add_ext(t2,t2,xQ,p);
                           add_ext(t3,xP4,xP4,p);
                           sub_ext(t2,t2,t3,p);
                           opposer(t2,t2,p);
                           mul_ext(xQres,t1,t2,p);
                           return;
                         }
  //calcul du coefficient A de la courbe de Montgomery
                         void get_A(comp_t A,comp_t xP,comp_t xQ,comp_t xR,mpz_t p)
                         {
                           comp_t t0,t1;
                           comp_init(t0);comp_init(t1);
                           add_ext(t1,xP,xQ,p);
                           mul_ext(t0,xP,xQ,p);
                           mul_ext(A,xR,t1,p);
                           add_ext(A,A,t0,p);
                           mul_ext(t0,t0,xR,p);
                           sub_ext_ui(A,A,1,p);
                           add_ext(t0,t0,t0,p);
                           add_ext(t1,t1,xR,p);
                           add_ext(t0,t0,t0,p);
                           mul_ext(A,A,A,p);
                           inv_ext(t0,t0,p);
                           mul_ext(A,A,t0,p);
                           sub_ext(A,A,t1,p);
                           return;
                         }

                       void get_yP_yQ_a_b(comp_t yP,comp_t yQ,comp_t a,comp_t b,comp_t xP, comp_t xQ,comp_t xR,
                                           mpz_t p)
                         {
                           comp_t t1,t2,xT,yT;
                           comp_init(t1);comp_init(t2);comp_init(xT);comp_init(yT);
                           get_A(a,xP,xQ,xR,p);
                           set_ext_ui(b,1);
                           mul_ext(t1,xP,xP,p);
                           mul_ext(t2,xP,t1,p);
                           mul_ext(t1,a,t1,p);
                           add_ext(t1,t2,t1,p);
                           add_ext(t1,t1,xP,p);
                           root_ext(yP,t1,p);
                           mul_ext(t1,xQ,xQ,p);
                         mul_ext(t2,xQ,t1,p);
                           mul_ext(t1,a,t1,p);
                           add_ext(t1,t2,t1,p);
                           add_ext(t1,t1,xQ,p);
                           root_ext(yQ,t1,p);
                           opposer(yQ,yQ,p);
                           xADD(xT,yT,xP,yP,xQ,yQ,a,b,p);
                           if(mpz_cmp(xT->x,xR->x)!=0||mpz_cmp(xT->y,xR->y)!=0)
                           opposer(yQ,yQ,p);
                           return;
                         }
                         void curve_3_iso(comp_t a_res,comp_t b_res,comp_t xP3,comp_t a,comp_t b,mpz_t p)
                          {
                            comp_t t1,t2;
                            comp_init(t1);comp_init(t2);
                            mul_ext(t1,xP3,xP3,p);
                            mul_ext(b_res,b,t1,p);
                            add_ext(t1,t1,t1,p);
                            add_ext(t2,t1,t1,p);
                            add_ext(t1,t1,t2,p);

                            set_ext_ui(t2,6);
                            sub_ext(t1,t1,t2,p);
                            mul_ext(t2,a,xP3,p);
                            sub_ext(t1,t2,t1,p);
                            mul_ext(a_res,t1,xP3,p);
                            return;
                         }

                        void point_3_iso(comp_t xQres,comp_t yQres,comp_t xQ,comp_t yQ,comp_t xP3,
                                         comp_t a, comp_t b,mpz_t p)
                            {
                              comp_t t1,t2,t3,t4;
                              comp_init(t1);comp_init(t2);comp_init(t3);comp_init(t4);
                              mul_ext(t1,xQ,xQ,p);
                              mul_ext(t1,t1,xP3,p);
                              mul_ext(t2,xP3,xP3,p);
                              mul_ext(t2,xQ,t2,p);
                              add_ext(t3,t2,t2,p);
                              add_ext(t2,t2,t3,p);
                              sub_ext(t1,t1,t2,p);
                              add_ext(t1,t1,xQ,p);
                              add_ext(t1,t1,xP3,p);
                              sub_ext(t2,xQ,xP3,p);
                              inv_ext(t2,t2,p);
                              mul_ext(t3,t2,t2,p);
                              mul_ext(t2,t2,t3,p);
                              mul_ext(t4,xQ,xP3,p);
                              sub_ext_ui(t4,t4,1,p);
                              mul_ext(t1,t4,t1,p);
                              mul_ext(t1,t1,t2,p);
                              mul_ext(t2,t4,t4,p);
                              mul_ext(t2,t2,t3,p);
                              mul_ext(xQres,xQ,t2,p);
                              mul_ext(yQres,yQ,t1,p);
                              return;
                            }

                            void iso_3_e(comp_t a_res,comp_t b_res,comp_t res_xP1,comp_t res_yP1,comp_t res_xP2,
                                        comp_t res_yP2,comp_t xS,comp_t yS,comp_t xP1,comp_t yP1,comp_t xP2,comp_t yP2,
                                         comp_t a,comp_t b,int e3,mpz_t p)
                                {
                                  int e=e3-1;
                                  int i=e;
                                  comp_t xT,yT;
                                  comp_init(xT);comp_init(yT);
                                  while(i>0){
                                         xTPLe(xT,yT,xS,yS,a,b,e,p);
                                         curve_3_iso(a_res,b_res,xT,a,b,p);
                                         point_3_iso(xS,yS,xS,yS,xT,a,b,p);
                                         i=i-1;
                                           }
                                           point_3_iso(res_xP1,res_yP1,xP1,yP1,xT,a,b,p);
                                           point_3_iso(res_xP2,res_yP2,xP2,yP2,xT,a,b,p);
                                     return;
                                   }
//calcul de la clé publique d'Alice
void isogen2(comp_t xP3_r,comp_t xQ3_r,comp_t xR3_r,mpz_t sk2,int e2,int e3,comp_t xP2,comp_t yP2,comp_t xQ2,comp_t yQ2,comp_t xP3,comp_t yP3,
     comp_t xQ3,comp_t yQ3,comp_t a,comp_t b,mpz_t p)
   {
      comp_t a_prime,b_prime;
      comp_t xS,yS,yP3_r,yQ3_r;
      comp_init(xS);comp_init(yS);comp_init(yP3_r);
      comp_init(yQ3_r);comp_init(a_prime);comp_init(b_prime);

      char m[500];
      mpz_get_str(m,2,sk2);
            double_and_add(xS,yS,m,xQ2,yQ2,a,b,p);
      xADD(xS,yS,xP2,yP2,xS,yS,a,b,p);
      iso_2_e(a_prime,b_prime,xP3_r,yP3_r,xQ3_r,yQ3_r,xS,yS,xP3,yP3,xQ3,yQ3,
                a,b,e2,p);

     get_xR(xR3_r,a,b,xP3_r,yP3_r,xQ3_r,yQ3_r,p);

      return;
  }
//calcul de la clé publique de Bob
void isogen3(comp_t xP2_r,comp_t xQ2_r,comp_t xR2_r,mpz_t sk3,int e2,int e3,comp_t xP2,comp_t yP2,comp_t xQ2,comp_t yQ2,comp_t xP3,comp_t yP3,
                           comp_t xQ3,comp_t yQ3,comp_t a,comp_t b,mpz_t p)
                         {
                            comp_t a_prime,b_prime;
                            comp_t xS,yS,yP2_r,yQ2_r;
                            comp_init(a_prime);comp_init(b_prime);
                            comp_init(yP2_r);comp_init(yQ2_r);
                            comp_init(xS);comp_init(yS);
                           char m[500];

                           mpz_get_str(m,2,sk3);
                          double_and_add(xS,yS,m,xQ3,yQ3,a,b,p);
                          xADD(xS,yS,xP3,yP3,xS,yS,a,b,p);
                          iso_3_e(a_prime,b_prime,xP2_r,yP2_r,xQ2_r,yQ2_r,xS,yS,xP2,yP2,xQ2,yQ2,
                            a,b,e3,p);

                          get_xR(xR2_r,a,b,xP2_r,yP2_r,xQ2_r,yQ2_r,p);
                          return;
                        }

   // calcul du secret d4Alice ,c'est le secret partagé
   void isoex2(comp_t j2,mpz_t sk2,comp_t xP_pk3,comp_t xQ_pk3,comp_t xR_pk3,int e2,int e3,
               comp_t a,comp_t b,mpz_t p)
       {
         comp_t op1,op2,op3,op4;
         comp_t  yP_pk3,yQ_pk3,xS,yS;
         comp_init(op1);comp_init(op2);comp_init(op3);comp_init(op4);
         comp_init(yP_pk3);comp_init(yQ_pk3);comp_init(xS);comp_init(yS);
         char m[500];
         mpz_get_str(m,2,sk2);
         set_ext_ui(op1,0);
         set_ext_ui(op2,0);
         set_ext_ui(op3,0);
         set_ext_ui(op4,0);
         get_yP_yQ_a_b(yP_pk3,yQ_pk3,a,b,xP_pk3,xQ_pk3,xR_pk3,p);
         double_and_add(xS,yS,m,xQ_pk3,yQ_pk3,a,b,p);
         xADD(xS,yS,xP_pk3,yP_pk3,xS,yS,a,b,p);
         iso_2_e(a,b,op1,op2,op3,op4,xS,yS,op1,op2,op3,op4,a,b,e2,p);
         j_inv(j2,a,b,p);

          }
//calcul du secret de Bob

void isoex3(comp_t j3,mpz_t sk3,comp_t xP_pk2,comp_t xQ_pk2,comp_t xR_pk2,int e2,int e3,
                   comp_t a,comp_t b,mpz_t p)
      {
      comp_t op1,op2,op3,op4;
      comp_t xS,yS,yP_pk2,yQ_pk2;
      comp_init(op1);comp_init(op2);comp_init(op3);comp_init(op4);
      comp_init(yP_pk2);comp_init(yQ_pk2);comp_init(xS);comp_init(yS);
      char m[500];
      mpz_get_str(m,2,sk3);
      set_ext_ui(op1,0);
      set_ext_ui(op2,0);
      set_ext_ui(op3,0);
      set_ext_ui(op4,0);
      get_yP_yQ_a_b(yP_pk2,xQ_pk2,a,b,xP_pk2,xQ_pk2,xR_pk2,p);
      double_and_add(xS,yS,m,xQ_pk2,yQ_pk2,a,b,p);
      xADD(xS,yS,xP_pk2,yP_pk2,xS,yS,a,b,p);
      iso_3_e(a,b,op1,op2,op3,op4,xS,yS,op1,op2,op3,op4,a,b,e3,p);
      j_inv(j3,a,b,p);
      return;
   }
