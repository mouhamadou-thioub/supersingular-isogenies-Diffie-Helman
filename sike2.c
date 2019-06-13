#include <math.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
//structure permettant de définir les éléments de Fp2
typedef struct{
    mpz_t x;
    mpz_t y;
   }comp;
typedef comp comp_t[1];
//les différentes procédures qui ont été définies plus bas
void comp_init(comp_t);
void add_ext(comp_t ,comp_t ,comp_t ,mpz_t );
void mul_ext(comp_t,comp_t,comp_t,mpz_t);
void mul_ext_ui(comp_t,comp_t,int,mpz_t );
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
void xDBL_bis(comp_t,comp_t,comp_t,comp_t,comp_t , comp_t ,mpz_t );
void xDBLe_bis(comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,int,mpz_t );
void xDBLADD(comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t ,comp_t,comp_t,comp_t,mpz_t p);
void xTPL_bis(comp_t,comp_t,comp_t ,comp_t ,comp_t , comp_t,mpz_t );
void xTPLe_bis(comp_t ,comp_t,comp_t,comp_t,comp_t,comp_t,int ,mpz_t );
void ladder3(comp_t,comp_t,char* ,comp_t,comp_t,comp_t ,comp_t,mpz_t );
void j_inv_bis(comp_t,comp_t,comp_t,mpz_t);
void get_A_bis(comp_t,comp_t,comp_t,comp_t ,mpz_t);
void curve_4_iso_bis(comp_t,comp_t,comp_t ,comp_t ,comp_t ,comp_t ,comp_t ,mpz_t);
void point_4_iso_bis(comp_t ,comp_t,comp_t,comp_t ,comp_t ,comp_t,comp_t,mpz_t);
void curve_3_iso_bis(comp_t,comp_t ,comp_t ,comp_t ,comp_t ,comp_t ,mpz_t );
void point_3_iso_bis(comp_t,comp_t,comp_t ,comp_t ,comp_t ,comp_t ,mpz_t );
void iso_2_e_bis(comp_t,comp_t,comp_t ,comp_t ,comp_t ,comp_t ,comp_t ,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,
                    comp_t,comp_t,comp_t,comp_t,int,mpz_t );

void iso_3_e_bis(comp_t,comp_t,comp_t ,comp_t ,comp_t ,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t,comp_t ,comp_t,
                  comp_t ,comp_t,comp_t,comp_t,int ,mpz_t );

void isogen2_bis(comp_t,comp_t,comp_t,mpz_t ,int ,int ,comp_t ,comp_t ,comp_t ,comp_t,comp_t,comp_t,mpz_t );
void isogen3_bis(comp_t,comp_t,comp_t,mpz_t ,int ,int ,comp_t,comp_t ,comp_t ,comp_t ,comp_t ,comp_t ,mpz_t );
void isoex2_bis(comp_t,comp_t,comp_t ,comp_t ,mpz_t ,int,mpz_t);
 void isoex3_bis(comp_t,comp_t,comp_t,comp_t,mpz_t ,int,mpz_t);
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
 comp_t xP2,xR2,xQ2,yQ2;//les coordonnés des points pubiques d'Alice
 comp_t xP3,xR3,xQ3,yQ3;//les coordonnées des points publiques de Bob
 mpz_init(xP2->x);mpz_init(xP2->y);
 mpz_init(xQ2->x);mpz_init(xQ2->y);
 mpz_init(xP3->x);mpz_init(xP3->y);
 mpz_init(xQ3->x);mpz_init(xQ3->y);
 mpz_init(xR2->x);mpz_init(xR2->y);mpz_init(xR3->x);mpz_init(xR3->y);
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

 printf("\nEntrer l'abscisse du point Q2 choisit par Alice:");
 gmp_scanf("%Zd %Zd",&x,&y);
 mpz_set(xQ2->x,x);mpz_set(xQ2->y,y);

 printf("\nEntrer l'abscisse  du point P3 choisit par Bob :");
 gmp_scanf("%Zd %Zd",&x,&y);
 mpz_set(xP3->x,x);mpz_set(xP3->y,y);


 printf("\nEntrer l'abscisse  du point Q3 choisit par Bob :");
 gmp_scanf("%Zd %Zd",&x,&y);
 mpz_set(xQ3->x,x);mpz_set(xQ3->y,y);

 printf("\nEntrer l'abscisse  du point R2 choisit par Bob :");
 gmp_scanf("%Zd %Zd",&x,&y);
 mpz_set(xR2->x,x);mpz_set(xR2->y,y);


 printf("\nEntrer l'abscisse  du point R3 choisit par Bob :");
 gmp_scanf("%Zd %Zd",&x,&y);
 mpz_set(xR3->x,x);mpz_set(xR3->y,y);




 mpz_clear(x);mpz_clear(y);
 mpz_t puiss2,puiss3,prod,var;
 mpz_init(puiss2);mpz_init(puiss3);mpz_init(prod);mpz_init(var);
 mpz_set_ui(var,2);
 mpz_pow_ui(puiss2,var,e2);
 mpz_set_ui(var,3);
 mpz_pow_ui(puiss3,var,e3);
 mpz_mul(prod,puiss2,puiss3);
 mpz_set_ui(var,1);
 mpz_sub(p,prod,var);

 printf("\n******************************************************************************\n");
 printf(" \nGénérations de la clé publique d'Alice   :\n");
 printf("\n******************************************************************************\n");
 mpz_urandomm(sk2,state,puiss2);

 isogen2_bis(xP3_prime,xQ3_prime,xR3_prime,sk2,e2,e3,xP2,xQ2,xR2,xP3,xQ3,xR3,p);

 gmp_printf ("\nla clé publique d'Alice est pk2=(%Zd + %Zdi/%Zd + %Zdi/%Zd + %Zdi \n)",
             xP3_prime->x,xP3_prime->y,xQ3_prime->x,xQ3_prime->y,xR3_prime->x,xR3_prime->y);


 printf("\n******************************************************************************");
 printf(" \nGénérations de la clé publique de Bob  :\n");
 printf("\n******************************************************************************");

 mpz_urandomm(sk3,state,puiss3);

 isogen3_bis(xP2_prime,xQ2_prime,xR2_prime,sk3,e2,e3,xP2,xQ2,xR2,xP3,xQ3,xR3,p);


 gmp_printf ("la clé publique d'Alice est pk3=(%Zd + %Zdi/%Zd + %Zdi/%Zd + %Zdi \n)",
             xP2_prime->x,xP2_prime->y,xQ2_prime->x,xQ2_prime->y,xR2_prime->x,xR2_prime->y);

 printf("\n******************************************************************************\n");
 printf(" \nGénérations de la clé privée d'Alice   :\n");
 printf("\n******************************************************************************\n");

 comp_t j2;//le j_invariant de la courbe secrète d'Alice
 comp_init(j2);
 isoex2_bis(j2,xP3_prime,xQ3_prime,xR3_prime,sk2,e2,p);

 gmp_printf ("la clé privée d'Alice est j2=%Zd + %Zdi \n",j2->x,j2->y);

 printf("\n******************************************************************************\n");
 printf(" \nGénérations de la clé privée de Bob   :\n");
 printf("\n******************************************************************************\n");

 comp_t j3;//le j_invariant de la courbe secrète d'Alice
 comp_init(j3);
 isoex3_bis(j3,xP2_prime,xQ2_prime,xR2_prime,sk3,e3,p);

 gmp_printf ("la clé privée d'Alice est j3=%Zd + i%Zd \n",j3->x,j3->y);

 return 0;
 }




void comp_init(comp_t C){
mpz_init(C->x);
mpz_init(C->y);
}

void add_ext(comp_t res, comp_t a,comp_t b,mpz_t p)
 {
   mpz_add(res->x,a->x,b->y);
   mpz_add(res->y,a->x,b->y);
   mpz_mod(res->x,res->x,p);
   mpz_mod(res->y,res->y,p);
}

void sub_ext(comp_t res,comp_t a,comp_t b,mpz_t p)
  {
    mpz_sub(res->x,a->x,b->x);
    mpz_sub(res->y,a->y,b->y);
    mpz_mod(res->x,res->x,p);
    mpz_mod(res->y,res->y,p);
  }

  void mul_ext(comp_t res,comp_t t1 ,comp_t t2,mpz_t p)
   { mpz_t inter1,inter2;
     mpz_init(inter1);mpz_init(inter2);
     mpz_mul(res->x,t1->x,t2->x);
     mpz_mul(inter1,t1->y,t2->y);
     mpz_sub(res->x,res->x,inter1);
     mpz_mul(res->y,t1->x,t2->y);
     mpz_mul(inter2,t2->x,t1->y);
     mpz_sub(res->y,res->y,inter2);
     mpz_mod(res->x,res->x,p);
     mpz_mod(res->y,res->y,p);
   }

void mul_ext_ui(comp_t res,comp_t t,int b,mpz_t p)
    {
      mpz_t a;
      mpz_init(a);
      mpz_set_ui(a,b);
      mpz_mul(res->x,t->x,a);
      mpz_set(res->y,t->y);
    }
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


      void add_ext_ui(comp_t res,comp_t a,int b,mpz_t p)
         {
           comp_t c;
           comp_init(c);
           mpz_set_ui(c->x,b);
           mpz_set_ui(c->y,0);
           add_ext(res,a,c,p);
         }
     void sub_ext_ui(comp_t res,comp_t a,int b,mpz_t p)
        {
          comp_t c;
          comp_init(c);
          mpz_set_ui(c->x,b);
          mpz_set_ui(c->y,0);
          sub_ext(res,a,c,p);
        }



        void inv_ext(comp_t res,comp_t t,mpz_t p)
         { mpz_t a,b,module,d,un;
           mpz_init(a);mpz_init(b);mpz_init(module);mpz_init(d);mpz_init(un);
           mpz_set_ui(un,1);
           mpz_mul(a,t->x,t->x);
           mpz_mul(b,t->y,t->y);
           mpz_add(module,a,b);
           mpz_gcd(d,module,p);
           if(mpz_cmp(d,un)==0)
           mpz_invert(module,module,p);
           mpz_mul(res->x,module,a);
           mpz_mul(res->y,module,b);
           mpz_mod(res->x,res->x,p);
           mpz_mod(res->y,res->y,p);
       }
       void racine(mpz_t res,mpz_t a,mpz_t p){
            mpz_t p1,exp;
            mpz_init(p1);mpz_init(exp);
            mpz_add_ui(p1,p,1);
            mpz_divexact_ui(exp,p1,4);
            mpz_powm(res,a,exp,p);
       }

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
           if(mpz_legendre(delta,p)!=1){
             mpz_sub(t2,a,delta);
             mpz_invert(deux,deux,p);
             mpz_mul(t2,t2,deux);
             if(mpz_legendre(t2,p)!=1){
               racine(res->x,t2,p);
               mpz_invert(t2,res->x,p);
               mpz_mul(prod,t2,deux);
               mpz_mul(res->y,b,prod);
             }
             else {
               mpz_set_ui(deux,2);
               mpz_add(t2,a,delta);
               mpz_invert(deux,deux,p);
               mpz_mul(t2,t2,deux);
               racine(res->x,t2,p);
               mpz_invert(t2,res->x,p);
               mpz_mul(prod,t2,deux);
               mpz_mul(res->y,b,prod);

             }
           }
         }

void xDBL_bis(comp_t x2P,comp_t z2P,comp_t xP,comp_t zP,comp_t A, comp_t C,mpz_t p)
  {
    comp_t t0,t1;
    comp_init(t0);comp_init(t1);
    sub_ext(t0,xP,zP,p);
    add_ext(t1,xP,zP,p);
    mul_ext(t0,t0,t0,p);
    mul_ext(t1,t1,t1,p);
    mul_ext(t1,t1,t1,p);
    mul_ext(z2P,C,t0,p);
    mul_ext(x2P,z2P,t1,p);
    sub_ext(t1,t1,t0,p);
    mul_ext(t0,A,t1,p);
    add_ext(z2P,z2P,t0,p);
    add_ext(z2P,z2P,t1,p);
    return;
  }

  void xDBLe_bis(comp_t x_res,comp_t z_res,comp_t xP,comp_t zP,comp_t A,comp_t C,int e,mpz_t p)
    {
      comp_t x1,z1;
      comp_init(x1);comp_init(z1);
      set_ext(x1,xP);
      set_ext(z1,zP);
      int i;
      for(i=0;i<=e;i++){
        xDBL_bis(x1,z1,x1,z1,A,C,p);
      }
      set_ext(x_res,x1);
      set_ext(z_res,z1);
      return;
    }
  void xDBLADD(comp_t x2P,comp_t z2P,comp_t xPQ,comp_t zPQ,comp_t xP,comp_t zP,comp_t xQ,comp_t zQ,
                      comp_t xQ_P,comp_t zQ_P,comp_t A,comp_t C,mpz_t p ){

      comp_t t0,t1,t2;
      comp_init(t0);comp_init(t1);comp_init(t2);
      add_ext(t0,xP,zP,p);
      sub_ext(t1,xP,zP,p);
      mul_ext(x2P,t0,t0,p);
      sub_ext(t2,xQ,zQ,p);
      add_ext(xPQ,xQ,zQ,p);
      mul_ext(t0,t0,t2,p);
      mul_ext(z2P,t1,t1,p);
      mul_ext(t1,t1,xPQ,p);
      sub_ext(t2,x2P,z2P,p);
      mul_ext(x2P,x2P,z2P,p);
      mul_ext(xPQ,A,t2,p);
      sub_ext(zPQ,t0,t1,p);
      add_ext(z2P,xPQ,z2P,p);
      add_ext(xPQ,t0,t1,p);
      mul_ext(z2P,z2P,t2,p);
      mul_ext(zPQ,zPQ,zPQ,p);
      mul_ext(zPQ,xQ_P,zPQ,p);
      mul_ext(xPQ,zQ_P,xPQ,p);
      return;
  }

  void xTPL_bis(comp_t x3P,comp_t z3P,comp_t xP,comp_t zP,comp_t A, comp_t C,mpz_t p){

      comp_t t0,t1,t2,t3,t4,t5,t6;
      comp_init(t0);comp_init(t1);comp_init(t2);comp_init(t3);
      comp_init(t4);comp_init(t5);comp_init(t6);
      sub_ext(t0,xP,zP,p);
      mul_ext(t2,t0,t0,p);
      add_ext(t1,xP,zP,p);
      mul_ext(t3,t1,t1,p);
      add_ext(t4,t1,t0,p);
      sub_ext(t0,t1,t0,p);
      mul_ext(t1,t4,t4,p);
      sub_ext(t1,t1,t3,p);
      sub_ext(t1,t1,t2,p);
      mul_ext(t5,t3,A,p);
      mul_ext(t3,t5,t3,p);
      mul_ext(t6,t2,C,p);
      mul_ext(t2,t2,t6,p);
      sub_ext(t3,t2,t3,p);
      sub_ext(t2,t5,t6,p);
      mul_ext(t1,t2,t1,p);
      add_ext(t2,t3,t1,p);
      mul_ext(t2,t2,t2,p);
      mul_ext(x3P,t2,t4,p);
      sub_ext(t2,t3,t1,p);
      mul_ext(t1,t1,t1,p);
      mul_ext(z3P,t1,t0,p);
      return;

  }
 void xTPLe_bis(comp_t x_res,comp_t z_res,comp_t xP,comp_t zP,comp_t A,comp_t C,int e,mpz_t p)
   {
     comp_t x1,z1;
     comp_init(x1);comp_init(z1);
     set_ext(x1,xP);
     set_ext(z1,zP);
     int i;
     for(i=0;i<=e;i++){
       xTPL_bis(x1,z1,x1,z1,A,C,p);
     }
     set_ext(x_res,x1);
     set_ext(z_res,z1);
     return;
   }
void ladder3(comp_t x1,comp_t z1,char* m,comp_t xP,comp_t xQ,comp_t xQ_P,comp_t A,mpz_t p)
    {
      comp_t x0,z0,x2,z2,C;
      comp_init(x0);comp_init(z0);comp_init(x2);comp_init(z2);comp_init(C);
      set_ext(x0,xQ); set_ext(x1,xP);set_ext(x2,xQ_P);
      set_ext_ui(z0,1);set_ext_ui(z1,1);set_ext_ui(z2,1);set_ext_ui(C,1);
      int i;
      int l=strlen(m);
      for(i=0;i<l;i++)
        {
          if(m[i]==1)
           xDBLADD(x0,z0,x1,z1,x0,z0,x1,z1,x2,z2,A,C,p);
           else
            xDBLADD(x0,z0,x2,z2,x0,z0,x2,z2,x1,z1,A,C,p);
        }
    }


void j_inv_bis(comp_t j,comp_t A,comp_t C,mpz_t p)
   {
     comp_t t0,t1;
     comp_init(t0);comp_init(t1);
     mul_ext(j,A,A,p);
     mul_ext(t1,C,C,p);
     add_ext(t0,t1,t1,p);
     sub_ext(t0,j,t0,p);
     sub_ext(t0,t0,t1,p);
     sub_ext(j,t0,t1,p);
     mul_ext(t1,t1,t1,p);
     mul_ext(j,j,t1,p);
     add_ext(t0,t0,t0,p);
     add_ext(t0,t0,t0,p);
     mul_ext(t1,t0,t0,p);
     mul_ext(t0,t0,t1,p);
     add_ext(t0,t0,t0,p);
     add_ext(t0,t0,t0,p);
     inv_ext(j,j,p);
     mul_ext(j,t0,j,p);

   }

  void get_A_bis(comp_t A,comp_t xP,comp_t xQ,comp_t xQ_P,mpz_t p)
    {
      comp_t t0,t1;
      comp_init(t0);comp_init(t1);
      add_ext(t1,xP,xQ,p);
      mul_ext(t0,xP,xQ,p);
      mul_ext(A,xQ_P,t1,p);
      add_ext(A,A,t0,p);
      mul_ext(t0,t0,xQ_P,p);
      sub_ext_ui(A,A,1,p);
      add_ext(t0,t0,t0,p);
      add_ext(t1,t1,xQ_P,p);
      add_ext(t0,t0,t0,p);
      mul_ext(A,A,A,p);
      inv_ext(t0,t0,p);
      mul_ext(A,A,t0,p);
      sub_ext(A,A,t1,p);
    }
  void curve_4_iso_bis(comp_t A,comp_t C,comp_t k1,comp_t k2,comp_t k3,comp_t xP4,comp_t zP4,mpz_t p)
      {
        sub_ext(k2,xP4,zP4,p);
        add_ext(k3,xP4,zP4,p);
        mul_ext(k1,zP4,zP4,p);
        add_ext(k1,k1,k1,p);
        mul_ext(C,k1,k1,p);
        add_ext(k1,k1,k1,p);
        mul_ext(A,xP4,xP4,p);
        add_ext(A,A,A,p);
        mul_ext(A,A,A,p);
      }
void point_4_iso_bis(comp_t xr,comp_t zr,comp_t k1,comp_t k2,comp_t k3,comp_t xQ,comp_t zQ,mpz_t p)
    {
      comp_t t0,t1;
      comp_init(t0);comp_init(t1);
      add_ext(t0,xQ,zQ,p);
      sub_ext(t1,xQ,zQ,p);
      mul_ext(xQ,t0,k2,p);
      mul_ext(zQ,t1,k3,p);
      mul_ext(t0,t1,t0,p);
      mul_ext(t0,t0,k1,p);
      add_ext(t1,xQ,zQ,p);
      sub_ext(zQ,xQ,zQ,p);
      mul_ext(t1,t1,t1,p);
      mul_ext(zQ,zQ,zQ,p);
      add_ext(zQ,t0,t0,p);
      sub_ext(t0,zQ,t0,p);
      mul_ext(xr,xQ,t1,p);
      mul_ext(zr,zQ,t0,p);
          }
void curve_3_iso_bis(comp_t A1,comp_t A2,comp_t k1,comp_t k2,comp_t xP3,comp_t zP3,mpz_t p)
     {
       comp_t t0,t1,t2,t3,t4;
       comp_init(t0);comp_init(t1);comp_init(t2);comp_init(t3);comp_init(t4);
       sub_ext(k1,xP3,zP3,p);
       mul_ext(t0,k1,k1,p);
       add_ext(k2,xP3,zP3,p);
       mul_ext(t1,k1,k1,p);
       add_ext(t2,t0,t1,p);
       add_ext(t3,k1,k2,p);
       mul_ext(t3,t3,t3,p);
       sub_ext(t3,t3,t2,p);
       add_ext(t2,t1,t3,p);
       add_ext(t3,t3,t0,p);
       add_ext(t4,t3,t0,p);
       add_ext(t4,t4,t4,p);
       add_ext(t4,t1,t4,p);
       mul_ext(A2,t2,t4,p);
       add_ext(t4,t1,t2,p);
       add_ext(t4,t4,t4,p);
       add_ext(t4,t0,t4,p);
       mul_ext(t4,t3,t4,p);
       sub_ext(t0,t4,A2,p);
       add_ext(A1,A2,t0,p);
     }
void point_3_iso_bis(comp_t x,comp_t z,comp_t k1,comp_t k2,comp_t xQ,comp_t zQ,mpz_t p)
    {
      comp_t t0,t1,t2;
      comp_init(t0);comp_init(t1);comp_init(t2);
      add_ext(t0,xQ,zQ,p);
      sub_ext(t1,xQ,zQ,p);
      mul_ext(t2,k1,t0,p);
      mul_ext(t1,k2,t1,p);
      add_ext(t2,t0,t1,p);
      sub_ext(t0,t1,t0,p);
      mul_ext(t2,t2,t2,p);
      mul_ext(t0,t0,t0,p);
      mul_ext(x,xQ,t2,p);
      mul_ext(z,zQ,t0,p);
    }
    void iso_2_e_bis(comp_t A,comp_t C,comp_t xP1_r,comp_t zP1_r,comp_t xP2_r,comp_t zP2_r,comp_t xP3_r,comp_t zP3_r,
             comp_t xS,comp_t zS,comp_t xP1,comp_t zP1,comp_t xP2,comp_t zP2,
                 comp_t xP3,comp_t zP3,comp_t a,comp_t c,int e2,mpz_t p)
        {
          int e=e2-2;
          int i=e;
          comp_t xT,zT;
          comp_init(xT);comp_init(zT);
          comp_t k1,k2,k3;
          comp_init(k1);comp_init(k2);comp_init(k3);
          while(i>0){
                 xDBLe_bis(xT,zT,xS,zS,a,c,e,p);
                 curve_4_iso_bis(A,C,k1,k2,k3,xT,zT,p);
                 point_4_iso_bis(xS,zS,k1,k2,k3,xS,zS,p);
                 i=i-2;
                   }

                point_4_iso_bis(xP1_r,zP1_r,k1,k2,k3,xP1,zP1,p);
                point_4_iso_bis(xP2_r,zP2_r,k1,k2,k3,xP2,zP2,p);
                point_4_iso_bis(xP3_r,zP3_r,k1,k2,k3,xP3,zP3,p);
        }

        void iso_3_e_bis(comp_t A,comp_t C,comp_t xP1_r,comp_t zP1_r,comp_t xP2_r,comp_t zP2_r,comp_t xP3_r,comp_t zP3_r,
                 comp_t xS,comp_t zS,comp_t xP1,comp_t zP1,comp_t xP2,comp_t zP2,
                     comp_t xP3,comp_t zP3,comp_t a,comp_t c,int e3,mpz_t p)
            {
              int e=e3-1;
              int i=e;
              comp_t xT,zT;
              comp_init(xT);comp_init(zT);
              comp_t k1,k2,k3;
              comp_init(k1);comp_init(k2);comp_init(k3);
              while(i>0){
                     xTPLe_bis(xT,zT,xS,zS,a,c,e,p);
                     curve_3_iso_bis(A,C,k1,k2,xT,zT,p);
                     point_3_iso_bis(xS,zS,k1,k2,xS,zS,p);
                     i=i-1;
                       }

                    point_3_iso_bis(xP1_r,zP1_r,k1,k2,xP1,zP1,p);
                    point_3_iso_bis(xP2_r,zP2_r,k1,k2,xP2,zP2,p);
                    point_3_iso_bis(xP3_r,zP3_r,k1,k2,xP3,zP3,p);
            }
// calcul de la clé publique d'Alice
    void isogen2_bis(comp_t x1,comp_t x2,comp_t x3,mpz_t sk2,int e2,int e3,comp_t xP2,comp_t xQ2,comp_t xR2,
                     comp_t xP3,comp_t xQ3,comp_t xR3,mpz_t p)

        {
             comp_t a0,c0,a,c,zP3,zQ3,zR3,xS,zS;
             comp_init(a0);comp_init(c0);comp_init(a);comp_init(c);
             comp_init(zP3);comp_init(zQ3);comp_init(zR3);comp_init(xS);comp_init(zS);
             set_ext_ui(a,0);set_ext_ui(c,1);
             set_ext_ui(a0,1);set_ext_ui(c0,2);
             set_ext_ui(zP3,1);set_ext_ui(zQ3,1);set_ext_ui(zR3,1);
             char m[500];
             mpz_get_str(m,2,sk2);
             ladder3(xS,zS,m,xP2,xQ2,xR2,a,p);
             iso_2_e_bis(a0,c0,xP3,zP3,xQ3,zQ3,xR3,zR3,a0,c0,xS,zS,xP3,zP3,xQ3,zQ3,xR3,zR3,e2,p);
             set_ext(x1,xP3);
             set_ext(x2,xQ3);
             set_ext(x3,xR3);

        }

//calcul de la clé publique de Bob
        void isogen3_bis(comp_t x1,comp_t x2,comp_t x3,mpz_t sk3,int e2,int e3,comp_t xP2,comp_t xQ2,comp_t xR2,
                         comp_t xP3,comp_t xQ3,comp_t xR3,mpz_t p)

            {
                 comp_t a0,c0,a,c,zP3,zQ3,zR3,xS,zS;
                 comp_init(a0);comp_init(c0);comp_init(a);comp_init(c);
                 comp_init(zP3);comp_init(zQ3);comp_init(zR3);comp_init(xS);comp_init(zS);
                 set_ext_ui(a,0);set_ext_ui(c,1);
                 set_ext_ui(a0,2);set_ext_ui(c0,-2);
                 set_ext_ui(zP3,1);set_ext_ui(zQ3,1);set_ext_ui(zR3,1);
                 char m[500];
                 mpz_get_str(m,2,sk3);
                 ladder3(xS,zS,m,xP2,xQ2,xR2,a,p);
                 iso_3_e_bis(a0,c0,xP3,zP3,xQ3,zQ3,xR3,zR3,a0,c0,xS,zS,xP3,zP3,xQ3,zQ3,xR3,zR3,e2,p);
                 set_ext(x1,xP2);
                 set_ext(x2,xQ3);
                 set_ext(x3,xR3);

            }
//calcul du secret d'Alice,c'est la clé commune
  void isoex2_bis(comp_t j2,comp_t x1,comp_t x2,comp_t x3,mpz_t sk2,int e2,mpz_t p)
       {
         comp_t a,c,a0,c0,xS,zS;
         comp_t xop1,zop1,xop2,zop2,xop3,zop3;
         comp_init(xop1),comp_init(xop2);comp_init(xop3);
         comp_init(zop1),comp_init(zop2);comp_init(zop3);
         comp_init(a);comp_init(c);comp_init(a0);comp_init(c0);
         comp_init(xS);comp_init(zS);
         get_A_bis(a,x1,x2,x3,p);
         set_ext_ui(c,1);
         char m[500];
         mpz_get_str(m,2,sk2);
         ladder3(xS,zS,m,x1,x2,x3,a,p);
         add_ext_ui(a0,a,2,p);
         set_ext_ui(c0,4);
         iso_2_e_bis(a0,c0,xop1,zop1,xop2,zop2,xop3,zop3,xS,zS,xop1,zop1,xop2,zop2,xop3,zop3,a0,c0,e2,p);
         mul_ext_ui(a0,a0,4,p);
         mul_ext_ui(c0,c0,2,p);
         sub_ext(a,a0,c0,p);
         j_inv_bis(j2,a,c0,p);

       }
//calcul du secret de Bob
       void isoex3_bis(comp_t j3,comp_t x1,comp_t x2,comp_t x3,mpz_t sk3,int e3,mpz_t p)
            {
              comp_t a,c,a0,c0,xS,zS;
              comp_t xop1,zop1,xop2,zop2,xop3,zop3;
              comp_init(xop1),comp_init(xop2);comp_init(xop3);
              comp_init(zop1),comp_init(zop2);comp_init(zop3);
              comp_init(a);comp_init(c);comp_init(a0);comp_init(c0);
              comp_init(xS);comp_init(zS);
              get_A_bis(a,x1,x2,x3,p);
              set_ext_ui(c,1);
              char m[500];
              mpz_get_str(m,2,sk3);
              ladder3(xS,zS,m,x1,x2,x3,a,p);
              add_ext_ui(a0,a,2,p);
              sub_ext_ui(c0,a0,2,p);
              iso_3_e_bis(a0,c0,xop1,zop1,xop2,zop2,xop3,zop3,xS,zS,xop1,zop1,xop2,zop2,xop3,zop3,a0,c0,e3,p);
              add_ext(a,a0,c0,p);
              sub_ext(c,a0,c0,p);
              mul_ext_ui(a,a,2,p);
              mul_ext_ui(c,c,2,p);
              j_inv_bis(j3,a,c0,p);

            }
