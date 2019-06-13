#include<stdlib.h>
#include<unistd.h>
#include<stdio.h>
int main(int argc,char ** argv)
{
  int choix;
  printf("                   MENU OPTION\n");
  printf("           ---------------------------------------------\n");
  printf("           : 1.SIDH version simple sans stratégié:\n");
  printf("           :      :\n");
  printf("           :           :\n");
  printf("           :            :\n");
  printf("                                          :\n");
  printf("           :  2. SIDH avec stratégie optimale                            :\n");
  printf("           :----------------------------------------------\n");

   printf("Faites votre choix:\n");
  scanf("%d",&choix);
  if (choix==1)
  execlp("./test1", "./test1",NULL);
  if(choix==2)
  execlp("./test2", "./test2",NULL);
  return 0;
  }
