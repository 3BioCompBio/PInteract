#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>
#include <float.h>

#define min(a,b) (a<=b ? a : b)
#define max(a,b) (a> b ? a : b)

#define NONEXIST 999.0
#define MEMMAX 100
#define HBMAX 20

//#define CATDmax 4.5
//#define CATangmax1 2.0

//#define PiPiDmax 5.0
//#define PiPiangmax1 3.0

//#define SPiDmax 6.0
//#define SPiangmax1 2.0

char ttamines[21],atbase[9][3];
char atcharg[5][5][4],ataro[4][9][4],atsulfur[2][1][4];
char rescharg[6],resaro[5],resulfur[3],resbase[10];
char residus[21];
int namines,nlong,frobs[20],tfrobs;
int atchmax[5],atarmax[4],atbamax[10],atsulfurmax[2];

float Bcoord[5][1][2][3],BcoordN[5][1][2][3],Bcoordtot[5][1][9][3];
float Acoord[5][1][2][3],AcoordN[5][1][2][3],Acoordtot[5][1][9][3];
float Ccoord[10][1][2][3],CcoordN[10][1][2][3],Ccoordtot[10][1][9][3];
float Dcoord[2][1][2][3],DcoordN[2][1][2][3],Dcoordtot[2][1][9][3];

float CATDmax=4.5,CATangmax1=2.0;
float PiPiDmax=5.0,PiPiangmax1=3.0;
float SPiDmax=6.0,SPiangmax1=2.0;

/*-------------------------------------------*/
char lett(char* res)
/*-------------------------------------------*/
{
    int i,j;
    char aa[2];

    if(!strcmp(res,"ALA"))aa[0]='A';
    else if(!strcmp(res,"ARG"))aa[0]='R';
    else if(!strcmp(res,"ASN"))aa[0]='N';
    else if(!strcmp(res,"ASP"))aa[0]='D';
    else if(!strcmp(res,"CYS"))aa[0]='C';
    else if(!strcmp(res,"GLN"))aa[0]='Q';
    else if(!strcmp(res,"GLU"))aa[0]='E';
    else if(!strcmp(res,"GLY"))aa[0]='G';
    else if(!strcmp(res,"HIS"))aa[0]='H';
    else if(!strcmp(res,"ILE"))aa[0]='I';
    else if(!strcmp(res,"LEU"))aa[0]='L';
    else if(!strcmp(res,"LYS"))aa[0]='K';
    else if(!strcmp(res,"MET"))aa[0]='M';
    else if(!strcmp(res,"PHE"))aa[0]='F';
    else if(!strcmp(res,"PRO"))aa[0]='P';
    else if(!strcmp(res,"SER"))aa[0]='S';
    else if(!strcmp(res,"THR"))aa[0]='T';
    else if(!strcmp(res,"TRP"))aa[0]='W';
    else if(!strcmp(res,"TYR"))aa[0]='Y';
    else if(!strcmp(res,"VAL"))aa[0]='V';
    else if(!strcmp(res,"UNK"))aa[0]='X';
    else if(!strcmp(res,"GLX"))aa[0]='Z';
    else if(!strcmp(res,"ASX"))aa[0]='B';

    else if(!strcmp(res,"HOH"))aa[0]='w';

    else
    {
      j=0;
      for(i=0;i<3;i++)
      {
          if(res[i]==' '||res[j]=='+')continue;
          res[j]=res[i];
          j++;
      }
      res[j]='\0';
      if(res[1]=='\0')aa[0]=res[0]+'a'-'A';
      else if(res[2]=='\0'&&res[0]=='D')aa[0]=res[1]+'a'-'A';
      else if(!strcmp(res,"G5'"))aa[0]='g';
      else
      {
          aa[0]='x';
      }
      aa[1]='\0';
    }

    return(*aa);
    
}
/*-------------------------------------------*/
float dist2(float* c1,
            float* c2)
/*-------------------------------------------*/
{
    int i;
    float d;

    for(i=0;i<3;i++)
        if(c1[i]==NONEXIST||c2[i]==NONEXIST)
        {
            d=NONEXIST;
            return(d);
        }

    d=0;
    for(i=0;i<3;i++)
        d+=(c1[i]-c2[i])*(c1[i]-c2[i]);
    return(d);
}
/*-------------------------------------------*/
float dist(float* c1,
           float* c2)
/*-------------------------------------------*/
{
    float d;
    d=(float)(sqrt((double)dist2(c1,c2)));
    return(d);
}
/*-------------------------------------------*/
float norme(float* c)
/*-------------------------------------------*/
{
    int i;
    float d;
    d=0.0;
    for(i=0;i<3;i++)d+=c[i]*c[i];
    d=(float)sqrt((double)d);
    if(d==0.0)
    {
        //  printf("Vanishing norm\n");
    }
    return(d);
}
/*-------------------------------------------*/
float angle(float* c1,
            float* c2,
            int flag)
/*-------------------------------------------*/
{
    int i;
    float d;
    d=0;
    for(i=0;i<3;i++)d+=c1[i]*c2[i];
    d/=(norme(c1)*norme(c2));
    d=(float)fabs(acos((double)d));
    if(d>M_PI/2.0+0.01)
    {
        if(flag)
        {
            // printf("Error in angle %f>pi/2=%f\n!!!!!!!!!!!!!!!!!!\n",d,M_PI/2);
        }
        d=M_PI/2.0;
    }
    return(d);
}
/*-------------------------------------------*/
bool anglemax(int nbase,
              int q,
              float dis,
              float *amax,
              float radcyl)
/*-------------------------------------------*/
{
    int j,k;
    float tt[3],ttN[3],a,d;
  
    d=(float)sqrt((double)dis);

    *amax=(-999.0);
    for(j=q*4;j<6+q*3;j++)
    {
        for(k=0;k<3;k++)
            tt[k]=(Bcoordtot[nbase][0][j][k]-Bcoord[nbase][0][0][k])*radcyl;
        for(k=0;k<3;k++)ttN[k]=BcoordN[nbase][0][0][k];
        for(k=0;k<3;k++)ttN[k]*=d;
        for(k=0;k<3;k++)ttN[k]+=tt[k];
        a=M_PI/2.0-angle(tt,ttN,true);
        if(a>*amax)*amax=a;
    }
    return(true);
}
/*-------------------------------------------*/
bool anglemaxaro(float dis,
                 float *amax,
                 float* cot,
                 float* co,
                 float* coN,
                 float radcyl)
/*-------------------------------------------*/
{
    int k;
    float tt[3],ttN[3],d;
  
    d=(float)sqrt((double)dis);
 
    *amax=(-999.0);
    for(k=0;k<3;k++)tt[k]=(cot[k]-co[k])*radcyl;
    for(k=0;k<3;k++)ttN[k]=coN[k];
    for(k=0;k<3;k++)ttN[k]*=d;
    for(k=0;k<3;k++)ttN[k]+=tt[k];
    *amax=M_PI/2.0-angle(tt,ttN,true);
   
    return(true);
}
/*-------------------------------------------*/
void milieu(float co[][3],
            float cofi[],
            int num,
            int* fl)
/*-------------------------------------------*/
{
    int i,j,cc;
    char aux[7];
  
    for(j=0;j<3;j++)cofi[j]=0.0;
    aux[6]='\0';
    cc=0;
    for(i=0;i<num;i++)
    {
        if(co[i][0]!=NONEXIST&&co[i][1]!=NONEXIST&&co[i][2]!=NONEXIST)
        {
            for(j=0;j<3;j++)cofi[j]+=co[i][j];
            cc++;
        }
    }
    if(cc!=0)
    {
        for(j=0;j<3;j++)cofi[j]/=(float)cc;
    }
    if(cc==0)*fl=false;
    else *fl=true;
}
/*-------------------------------------------*/
bool normale(float comil[],
             float co[][3],
             int num,
             float coN[])
/*-------------------------------------------*/
{
    int i,j,k,cc;
    float cici[(6*5)/2][3],tt[3],norm,x1,x2,y1,y2,z1,z2,d1,d2;
 
    if(num>6)
    {
        printf("ERROR 1\n");
        // printf("Strange, num=%d>6\n",num);
        exit(1);
    }
    cc=0;
    for(i=0;i<num;i++)
    {
        if(co[i][0]==NONEXIST||co[i][1]==NONEXIST||co[i][2]==NONEXIST)
            continue;
        cc++;
    }
    if(cc<2)
    {
        //  printf("Impossible to calculate the norm, I have only 2 points\n");
        for(k=0;k<3;k++)coN[k]=NONEXIST;
        return(false);
    }
    cc=0;
    for(i=0;i<num;i++)
    {
        if(co[i][0]==NONEXIST||co[i][1]==NONEXIST||co[i][2]==NONEXIST)
            continue;
        x1=co[i][0]-comil[0];
        y1=co[i][1]-comil[1];
        z1=co[i][2]-comil[2];
        for(j=i+1;j<num;j++)
        {
            if(co[j][0]==NONEXIST||co[j][1]==NONEXIST||co[j][2]==NONEXIST)
                continue;
            x2=co[j][0]-comil[0];
            y2=co[j][1]-comil[1];
            z2=co[j][2]-comil[2];
            cici[cc][0]=y1*z2-y2*z1;
            cici[cc][1]=z1*x2-z2*x1;
            cici[cc][2]=x1*y2-x2*y1;
            if((norm=norme(cici[cc]))==0)
                continue;
            for(k=0;k<3;k++)cici[cc][k]/=norm;
            cc++;
        }
    }
    if(cc==0)
    {
        //   printf("Impossible to calculate the norm, vanishing\n");
        for(k=0;k<3;k++)coN[k]=NONEXIST;
        return(false);
        }
    for(i=1;i<cc;i++)
    {
        d1=dist2(cici[i],cici[0]);
        for(j=0;j<3;j++)tt[j]=(-cici[i][j]);
        d2=dist2(tt,cici[0]);
        if(d2<d1)
            for(j=0;j<3;j++)cici[i][j]=tt[j];
    }
    for(j=0;j<3;j++)coN[j]=0.0;
    for(i=0;i<cc;i++)
        for(j=0;j<3;j++)
            coN[j]+=cici[i][j];
    for(j=0;j<3;j++)coN[j]/=(float)cc;
    norm=norme(coN);
    for(k=0;k<3;k++)coN[k]/=norm;
    return(true);
}
/*-------------------------------------------*/
void read_protADN_length(char* fichin,
                  int* nres)
/*---------------------------*/
{
    char no[20],noold[20],nomres[10],cr[2];
    char mot[256],ligne[2560];
    int jj;
    FILE *seqin;
    
    
    seqin=fopen(fichin,"r");
    if(seqin==NULL)
    {
        printf("*** '%s' not existing ***\n",fichin);
        printf("ERROR 2\n");
        exit(2);
    }
    
    mot[0]='\0';
    strcpy(noold,"vide");
    strcpy(no,"xxxxxx");
    nomres[3]='\0';
    *nres=0;
    
    while(fgets(ligne,2560,seqin)!=NULL)
    {
        sscanf(ligne,"%s",mot);
        strncpy(no,&ligne[21],6);
        if(!strncmp(ligne,"ENDMDL",6))break;
        strncpy(nomres,&ligne[17],3);
        cr[0]=lett(nomres);
        if(strchr(residus,cr[0])==NULL)continue;
        if((!strcmp(mot,"ATOM")||!strncmp(mot,"HETATM",6))&&strncmp(nomres,"HOH",3)&&strncmp(no,noold,6))
        {
            if(cr[0]=='x')
                (*nres)+=5;
            else
                (*nres)++;
            strcpy(noold,no);
        }
    }
    
    (*nres)++;

    fclose(seqin);
    
}
/*-------------------------------------------*/
void read_protADN(char* fichin,
                  int* nres,
                  char* seq,
                  char* longseq,
                  char* prefixes,
                  char* hist,
                  float Icoord[][2][3],
                  float IcoordN[][2][3],
                  float Icoordtot[][9][3],
                  char* flagADN,
                  char* name)
/*---------------------------*/
{
    int i,j,k,ok,jj,ii,kk,ll,bon,fl,mm,nn,ibase,lat,i1,i2,i3,i4,bb;
    float x,y,z,co[9][3],cotmp[3];
    char st[256],at[20],no[20],noold[20],hdeb[7],hfin[7];
    char cr[2],prefix[2],prefixold[2];
    char mot[256],ligne[2560],nomres[10],aux[7],aux1[4],aux2[5];
    FILE *seqin;
    
    hdeb[0]='\0';
    hfin[0]='\0';
    aux[6]='\0';
    aux1[3]='\0';
    aux2[4]='\0';
    no[6]='\0';
    name[0]='\0';
    jj=0;
    
    mot[0]='\0';
    strcpy(noold,"vide");
    strcpy(no,"xxxxxx");
    
    seqin=fopen(fichin,"r");
    if(seqin==NULL)
    {
        printf("*** '%s' not existing ***\n",fichin);
        printf("ERROR 3\n");
        exit(3);
    }
    
    while(fgets(ligne,2560,seqin)!=NULL)
    {
        sscanf(ligne,"%s",mot);
        if(!strcmp(mot,"TITLE"))strncat(name,&ligne[10],60);
        if(!strcmp(mot,"COMPND")&&name[0]=='\0')strncat(name,&ligne[10],60);
        if(!strcmp(mot,"ATOM")||!strcmp(mot,"HETATM"))break;
    }
    if(strcmp(mot,"ATOM")&&strcmp(mot,"HETATM"))
    {
        printf("*** '%s' incomplete ***\n",fichin);
        //  printf("**** '%s' ****\n",mot);
        printf("ERROR 4\n");
        exit(4);
    }
    
    ok=false;
    bon=false;
    cr[0]='\0';
    strcpy(no,"xxxxxx");
    strcpy(noold,"vide");
    nomres[3]='\0';
    prefix[0]='\0';
    prefixold[0]='\0';
    ibase=(-1);
    
    if(name[0]=='\0')strcpy(name,"NoName");
    else
    {
        j=(int)strlen(name);
        for(i=1;i<j;i++)
        {
            if(name[i]==' '&&name[i-1]==' ')
            {
                for(k=i;k<j;k++)name[k]=name[k+1];
                j--;i--;
                continue;
            }
        }
        if(name[j-1]==' ')name[j-1]='\0';
    }

    while(true)
    {
        st[0]='\0';
        at[0]='\0';
        sscanf(ligne,"%s",st);
        if(!strncmp(st,"ENDMDL",6))break;
        if((!strncmp(st,"ATOM",4)||!strncmp(st,"HETATM",6))&&strncmp(&ligne[17],"HOH",3))
        {
            strncpy(aux2,&ligne[12],4);
            sscanf(aux2,"%s",at);
            strncpy(nomres,&ligne[17],3);
            if(ligne[16]!=' '&&ligne[16]!='A'&&ligne[16]!='1')
            {
                if(fgets(ligne,2560,seqin)==NULL)break;
                continue;
            }
            strncpy(no,&ligne[21],6);
            if(lett(nomres)=='x')
            {
                lat=(int)strlen(at);
                i1=(int)(strstr(at,"C5A")-at);
                i2=(int)(strstr(at,"C5M")-at);
                i3=(int)(strstr(at,"CN")-at);
                if((i1>=0&&i1<lat&&prefixold[0]!='A')||
                   (i2>=0&&i2<lat&&prefixold[0]!='M')||
                   (i3>=0&&i3<lat&&prefixold[0]!='C'&&prefixold[0]!='N'))
                    prefix[0]=prefixold[0];
                else if(at[1]>='0'&&at[1]<='9')
                    prefix[0]=at[2];
                else if(at[2]>='0'&&at[2]<='9')
                    prefix[0]=at[0];
                else prefix[0]='\0';
                if(prefix[0]=='*')prefix[0]='\0';
            }
            else prefix[0]='\0';
            sscanf(&ligne[30],"%f %f %f",&x,&y,&z);
            if(!ok)
            {
                if(hdeb[0]=='\0')strcpy(hdeb,no);
                if(strcmp(no,hdeb))
                {
                    if(fgets(ligne,2560,seqin)==NULL)break;
                    continue;
                }
                ok=true;
                jj=(-1);
            }
            if(!strcmp(noold,no)&&!strncmp(prefixold,prefix,1)&&cr[0]=='x')
                if(!strncmp(&longseq[jj*3],"CA",2))
                {
                    strncpy(&longseq[jj*3],nomres,3);
                    if(longseq[jj*3+1]=='\0')longseq[jj*3+1]=' ';
                    if(longseq[jj*3+2]=='\0')longseq[jj*3+2]=' ';
                    longseq[jj*3+3]='\0';
                }
            
            if(strcmp(noold,no)||strncmp(prefixold,prefix,1))
            {
                if(bon&&cr[0]!='\0')
                {
                    if(cr[0]=='x')
                    {
                        for(i=0;i<9;i++)
                        {
                            for(j=0;j<3;j++)if(co[i][j]==NONEXIST)break;
                            if(j<3)break;
                        }
                        strncpy(aux1,&longseq[jj*3],3);
                        if(i<6||i==7||i==8)
                        {
                            jj--;
                            bon=false;
                        }
                        else
                        {
                            if(ibase<0)
                            {
                                //   printf("error, ibase=%d\n",ibase);
                                //   printf("### %s  %s\n",fichin,aux1);
                                if(i==9)
                                {
                                    seq[jj]='y';
                                    cr[0]='y';
                                }
                                else
                                {
                                    for(bb=0;bb<strlen(fichin)-3;bb++)
                                        if(!strncmp(&fichin[bb],"1b6t",4))break;
                                    if(bb<strlen(fichin)-3)
                                    {
                                        seq[jj]='a';
                                        cr[0]='a';
                                        ibase=(int)(strchr(resbase,'a')-resbase);
                                        // printf("### Corrected\n");
                                    }
                                }
                            }
                            else
                            {
                                seq[jj]=resbase[ibase];
                                cr[0]=resbase[ibase];
                            }
                        }
                    }
                    if(bon)
                    {
                        if(cr[0]=='R')ii=2;
                        else if(cr[0]=='Q'||cr[0]=='N')ii=1;
                        else if(cr[0]=='K')ii=1;
                        else if(cr[0]=='H')ii=5;
                        else if(cr[0]=='Y'||cr[0]=='W'||cr[0]=='F')ii=6;
                        else if(cr[0]=='M'||cr[0]=='C')ii=1;
                        else if(cr[0]=='c'||cr[0]=='t'||cr[0]=='g'||cr[0]=='a'
                                ||cr[0]=='u'||cr[0]=='x'||cr[0]=='y')
                            ii=6;
                        else
                        {
                          //  printf("What is that value of cr %c\n",cr[0]);
                            printf("ERROR 5\n");
                            exit(5);
                        }
                        for(kk=0;kk<2;kk++)
                            for(ll=0;ll<3;ll++)
                            {
                                Icoord[jj][kk][ll]=0.0;
                                IcoordN[jj][kk][ll]=0.0;
                            }
                        for(kk=0;kk<9;kk++)
                            for(ll=0;ll<3;ll++)
                                Icoordtot[jj][kk][ll]=co[kk][ll];
                        milieu(co,Icoord[jj][0],ii,&fl);
                        if(fl==0)jj--;
                        else
                        {
                            if(cr[0]=='c'||cr[0]=='t'||cr[0]=='g'||cr[0]=='a'||
                               cr[0]=='u'||cr[0]=='x'||cr[0]=='y'||
                               cr[0]=='Y'||cr[0]=='W'||cr[0]=='F'||cr[0]=='H')
                                normale(Icoord[jj][0],co,ii,IcoordN[jj][0]);
                            else if(cr[0]=='Q'||cr[0]=='N')
                                normale(Icoord[jj][0],&co[1],2,IcoordN[jj][0]);
                            else if(cr[0]=='R')
                                normale(Icoord[jj][0],&co[1],2,IcoordN[jj][0]);
                            if(cr[0]=='g'||cr[0]=='a'||cr[0]=='y'||cr[0]=='W')
                            {
                                milieu(&co[4],Icoord[jj][1],5,&fl);
                                normale(Icoord[jj][1],&co[4],5,cotmp);
                                for(ll=0;ll<3;ll++)IcoordN[jj][1][ll]=(-cotmp[ll]);
                            }
                        }
                        bon=false;
                    }
                }
                if(!strcmp(noold,hfin))
                {
                    strcpy(no,noold);
                    break;
                }
                strcpy(noold,no);
                prefixold[0]=prefix[0];
                
                cr[0]=lett(nomres);
                if(cr[0]=='x')flagADN[jj+1]='2';
                else if(cr[0]>='a'&&cr[0]<='z')flagADN[jj+1]='1';
                else flagADN[jj+1]='0';
                if(strchr(residus,cr[0])!=NULL)bon=true;
                if(bon)
                {
                    jj++;
                    seq[jj]=cr[0];
                    strncpy(&hist[jj*6],no,6);
                    strncpy(&longseq[jj*3],nomres,3);
                    if(longseq[jj*3+1]=='\0')longseq[jj*3+1]=' ';
                    if(longseq[jj*3+2]=='\0')longseq[jj*3+2]=' ';
                    longseq[jj*3+3]='\0';
                    prefixes[jj]=prefix[0];
                    for(i=0;i<9;i++)
                        for(j=0;j<3;j++)co[i][j]=NONEXIST;
                    ibase=(-1);
                }
            }
            if(bon&&flagADN[jj]=='0')
            {
                ii=(-1);
                for(mm=0;mm<strlen(resaro);mm++)
                    if(cr[0]==resaro[mm])break;
                if(mm<strlen(resaro))
                {
                    for(nn=0;nn<atarmax[mm];nn++)
                        if(!strcmp(ataro[mm][nn],at))break;
                    if(nn<atarmax[mm])ii=nn;
                }
                else
                {
                    for(mm=0;mm<strlen(rescharg);mm++)
                        if(cr[0]==rescharg[mm])break;
                    if(mm<strlen(rescharg))
                    {
                        for(nn=0;nn<atchmax[mm];nn++)
                            if(!strcmp(atcharg[mm][nn],at))break;
                        if(nn<atchmax[mm])ii=nn;
                    }
                    else
                    {
                        for(mm=0;mm<strlen(resulfur);mm++)
                            if(cr[0]==resulfur[mm])break;
                        if(mm<strlen(resulfur))
                        {
                            for(nn=0;nn<atsulfurmax[mm];nn++)
                                if(!strcmp(atsulfur[mm][nn],at))break;
                            if(nn<atsulfurmax[mm])ii=nn;
                        }
                    }
                }
                if(ii>=0)
                {
                    if(co[ii][0]!=NONEXIST||co[ii][1]!=NONEXIST||co[ii][2]!=NONEXIST)
                    {
                      //  strncpy(aux,&hist[jj*6],6);
                      //  printf(
                      //         "**** Strange, 2 Icoordinates for atom %d of residue %s in %s ****\n",
                      //         i,aux,fichin);
                    }
                    co[ii][0]=x;
                    co[ii][1]=y;
                    co[ii][2]=z;
                }
            }
            else if(bon&&flagADN[jj]=='1')
            {
                for(ii=0;ii<9;ii++)
                    if(!strcmp(at,atbase[ii]))break;
                if(ii==9)ii=(-1);
                if(ii>=0)
                {
                    if(co[ii][0]!=NONEXIST||co[ii][1]!=NONEXIST||co[ii][2]!=NONEXIST)
                    {
                      // strncpy(aux,&hist[jj*6],6);
                      //  printf(
                      //         "**** Strange, 2 coordinates for atom %d of DNA %s in %s ****\n",
                      //         i,aux,fichin);
                    }
                    co[ii][0]=x;
                    co[ii][1]=y;
                    co[ii][2]=z;
                }
            }
            else if(bon&&flagADN[jj]=='2')
            {
                lat=(int)strlen(at);
                i3=(int)(strchr(at,'*')-at);
                for(ii=0;ii<9;ii++)
                    if(!strncmp(at,atbase[ii],2)||!strncmp(&at[1],atbase[ii],2))
                        break;
                i4=(int)(strstr(at,"CN")-at);
                if(i3>=0&&i3<lat)ii=10;
                else if(i4>=0&&i4<lat&&prefix[0]!='C'&&prefix[0]!='N')ii=10;
                else
                {
                    i1=(int)(strstr(at,"C5A")-at);
                    i2=(int)(strstr(at,"C5M")-at);
                    if((i1>=0&&i1<lat&&prefix[0]!='A')||
                       (i2>=0&&i2<lat&&prefix[0]!='M'))
                    {
                        ibase=(int)(strchr(resbase,'t')-resbase);
                        ii=10;
                    }
                }
                if(ii==9)
                {
                    if(!strncmp(at,"N2",2)||!strncmp(&at[1],"N2",2))
                        ibase=(int)(strchr(resbase,'g')-resbase);
                    else if(!strncmp(at,"N6",2)||!strncmp(&at[1],"N6",2))
                        ibase=(int)(strchr(resbase,'a')-resbase);
                    else if(!strncmp(at,"N4",2)||!strncmp(&at[1],"N4",2))
                        ibase=(int)(strchr(resbase,'c')-resbase);
                    else if((!strncmp(at,"O4",2)||!strncmp(&at[1],"O4",2))&&ibase<0)
                        ibase=(int)(strchr(resbase,'u')-resbase);
                    else if(!strncmp(at,"C7",2)||!strncmp(&at[1],"C7",2))
                        ibase=(int)(strchr(resbase,'t')-resbase);
                }
                if(ii>=0&&ii<9)
                {
                    if(co[ii][0]!=NONEXIST||co[ii][1]!=NONEXIST||co[ii][2]!=NONEXIST)
                    {
                     //   strncpy(aux,&hist[jj*6],6);
                     //   printf("**** Strange, 2 Icoordinates for atom %d of DNA %s in %s ****\n",
                     //          i,aux,fichin);
                    }
                    co[ii][0]=x;
                    co[ii][1]=y;
                    co[ii][2]=z;
                }
            }
        }
        if(fgets(ligne,2560,seqin)==NULL)break;
    }
        
    if(bon)
    {
        if(cr[0]=='x')
        {
            for(i=0;i<9;i++)
            {
                for(j=0;j<3;j++)if(co[i][j]==NONEXIST)break;
                if(j<3)break;
            }
            strncpy(aux1,&longseq[jj*3],3);
            if(i<6||i==7||i==8)
            {
                jj--;
                bon=false;
            }
            else
            {
                if(ibase<0)
                {
                 //   printf("error, ibase=%d\n",ibase);
                 //   printf("### %s  %s\n",fichin,aux1);
                    if(i==9)
                    {
                        seq[jj]='y';
                        cr[0]='y';
                    }
                    else
                    {
                        for(bb=0;bb<strlen(fichin)-3;bb++)
                            if(!strncmp(&fichin[bb],"1b6t",4))break;
                        if(bb<strlen(fichin)-3)
                        {
                            seq[jj]='a';
                            cr[0]='a';
                            ibase=(int)(strchr(resbase,'a')-resbase);
                            //    printf("### Corrected\n");
                        }
                    }
                }
                else
                {
                    seq[jj]=resbase[ibase];
                    cr[0]=resbase[ibase];
                }
            }
        }
        if(bon)
        {
            if(cr[0]=='R')ii=2;
            else if(cr[0]=='Q'||cr[0]=='N')ii=1;
            else if(cr[0]=='K')ii=1;
            else if(cr[0]=='H')ii=5;
            else if(cr[0]=='Y'||cr[0]=='W'||cr[0]=='F')ii=6;
            else if(cr[0]=='C'||cr[0]=='M')ii=1;
            else if(cr[0]=='c'||cr[0]=='t'||cr[0]=='g'||cr[0]=='a'
                    ||cr[0]=='u'||cr[0]=='x'||cr[0]=='y')ii=6;
            else
            {
               // printf("What is that cr value %c\n",cr[0]);
                printf("ERROR 6\n");
                exit(6);
            }
            for(kk=0;kk<9;kk++)
                for(ll=0;ll<3;ll++)
                    Icoordtot[jj][kk][ll]=co[kk][ll];
            milieu(co,Icoord[jj][0],ii,&fl);
            if(fl==0)jj--;
            else
            {
                if(cr[0]=='c'||cr[0]=='t'||cr[0]=='g'||cr[0]=='a'||
                    cr[0]=='u'||cr[0]=='x'||cr[0]=='y'||
                    cr[0]=='Y'||cr[0]=='W'||cr[0]=='F'||cr[0]=='H')
                    normale(Icoord[jj][0],co,ii,IcoordN[jj][0]);
                else if(cr[0]=='Q'||cr[0]=='N')
                    normale(Icoord[jj][0],&co[1],2,IcoordN[jj][0]);
                else if(cr[0]=='R')
                    normale(Icoord[jj][0],&co[1],2,IcoordN[jj][0]);
                if(cr[0]=='g'||cr[0]=='a'||cr[0]=='y'||cr[0]=='W')
                {
                    milieu(&co[4],Icoord[jj][1],5,&fl);
                    normale(Icoord[jj][1],&co[4],5,cotmp);
                    for(ll=0;ll<3;ll++)IcoordN[jj][1][ll]=(-cotmp[ll]);
                }
            }
        }
    }
    if(hfin[0]=='\0')strcpy(hfin,no);
    if(strcmp(no,hfin))
    {
     //   printf("**** %s or %s not found in %s ****\n",hdeb,hfin,fichin);
     //   printf("**** Program unsuccessful ****\n");
        printf("ERROR 7\n");
        exit(7);
    }
    *nres=jj+1;
    seq[jj+1]='\0';
    longseq[(jj+1)*3]='\0';
    hist[(jj+1)*6]='\0';
    flagADN[jj+1]='\0';
        
    /*
        j=0;
        for(i=0;i<(*nres);i++)
        if(flagADN[i]=='1')j++;
        if(j<=3)
        {for(i=0;i<(*nres);i++)
        if(flagADN[i]=='1')flagADN[i]='2';
        }
    */
    
    fclose(seqin);
}
/*-------------------------------------------*/
void read_hbonds(char fich[],
                 char repert[],
                 int nres,
                 char seq[],
                 char hist[],
                 char hb_seq[][HBMAX],
                 char hb_hist[][6*HBMAX],
                 char hb_at[][7*HBMAX],
                 char flhbplus[])
/*-------------------------------------------*/
{
    int i,i1,i2,KK;
    char fichin[128],ligne[2560],aux1[7],aux2[7];
    char res1[4],res2[4],at1[4],at2[4];
    FILE *seqin;

    aux1[6]='\0';
    aux2[6]='\0';

    sprintf(fichin,"%s%s",repert,fich);
    i1=(int)(strrchr(fichin,'.')-fichin);
    if(i1<0||i1>=strlen(fichin))
    {
      //  printf("No . in %s\n",fichin);
        printf("ERROR 8\n");
        exit(8);
    }
    strcpy(&fichin[i1+1],"hb2");
    seqin=fopen(fichin,"r");
    flhbplus[0]=(char)true;
    if(seqin==NULL)
    {
        printf("Missing hbplus file %s - No H-bonds available\n",fichin);
        flhbplus[0]=(char)false;
        return;
    }

    for(i1=0;i1<8;i1++)
        if(fgets(ligne,2560,seqin)==NULL)
        {
            printf("Hbplus file %s incomplete - H-bonds ignored\n",fichin);
            flhbplus[0]=(char)false;
            return;
        }
    while(fgets(ligne,2560,seqin)!=NULL)
    {
        strncpy(aux1,ligne,6);
        strncpy(aux2,&ligne[14],6);
        if(aux1[0]=='-')aux1[0]=' ';
        if(aux2[0]=='-')aux2[0]=' ';
        if(aux1[5]=='-')aux1[5]=' ';
        if(aux2[5]=='-')aux2[5]=' ';
        if(aux1[1]=='-')
        {
            aux1[1]=' ';
            i1=false;
        }
        else i1=true;
        if(aux2[1]=='-')
        {
            aux2[1]=' ';
            i2=false;
        }
        else i2=true;
        for(i=1;i<5;i++)
        {
            if(aux1[i]>='1'&&aux1[i]<='9')
            {
                if(!i1)aux1[i-1]='-';
                break;
            }
            if(aux1[i]=='0')aux1[i]=' ';
        }
        for(i=1;i<5;i++)
        {
            if(aux2[i]>='1'&&aux2[i]<='9')
            {
                if(!i2)aux2[i-1]='-';
                break;
            }
            if(aux2[i]=='0')aux2[i]=' ';
        }
    
        for(i1=0;i1<nres;i1++)
            if(!strncmp(aux1,&hist[i1*6],6))break;
        for(i2=0;i2<nres;i2++)
            if(!strncmp(aux2,&hist[i2*6],6))break;
        for(i=0;i<4;i++)
        {
            res1[i]='\0';
            res2[i]='\0';
        }
        strncpy(res1,&ligne[6],3);
        strncpy(res2,&ligne[20],3);
        res1[0]=lett(res1);
        res2[0]=lett(res2);
        strncpy(at1,&ligne[10],3);
        strncpy(at2,&ligne[24],3);
        if(res1[0]=='w'||res2[0]=='w'||res1[0]=='X'||res2[0]=='X')continue;
        if(res1[0]=='o'||res2[0]=='o')continue;
        if(res1[0]>='A'&&res1[0]<='Z'&&res2[0]>='A'&&res2[0]<='Z')continue;
        if(res1[0]>='a'&&res1[0]<='z'&&res2[0]>='a'&&res2[0]<='z')continue;
        if(i1<nres)
        {
            KK=true;
            if(res1[0]=='x'&&(seq[i1]=='a'||seq[i1]=='c'||seq[i1]=='g'
                              ||seq[i1]=='t'||seq[i1]=='u'||seq[i1]=='x'||seq[i1]=='y'))
            {}
            else if(res1[0]!=seq[i1])
            {
              //  printf("%s : %c != %c\n",aux1,res1[0],seq[i1]);
                KK=false;
            }
            for(i=0;i<strlen(hb_seq[i1]);i++)
                if(!strncmp(&hb_hist[i1][i*6],aux2,6))break;
            if(i==strlen(hb_seq[i1])&&KK)
            {
                strncat(hb_seq[i1],&seq[i2],1);
                strncat(hb_hist[i1],aux2,6);
                strncat(hb_at[i1],at1,3);
                strcat(hb_at[i1],"-");
                strncat(hb_at[i1],at2,3);
                if(strlen(hb_seq[i1])>=HBMAX)
                {
                   // printf("Increase HBMAX\n");
                    printf("ERROR 9\n");
                    exit(9);
                }
            }
        }
        if(i2<nres)
        {
            KK=true;
            if(res2[0]=='x'&&(seq[i2]=='a'||seq[i2]=='c'||seq[i2]=='g'
                              ||seq[i2]=='t'||seq[i2]=='u'||seq[i2]=='x'||seq[i2]=='y'))
            {}
            else if(res2[0]!=seq[i2])
            {
              //  printf("%s : %c != %c\n",aux2,res2[0],seq[i2]);
                KK=false;
            }
            for(i=0;i<strlen(hb_seq[i2]);i++)
                if(!strncmp(&hb_hist[i2][i*6],aux1,6))break;
            if(i==strlen(hb_seq[i2])&&KK)
            {
                strncat(hb_seq[i2],&seq[i1],1);
                strncat(hb_hist[i2],aux1,6);
                strncat(hb_at[i2],at2,3);
                strcat(hb_at[i2],"-");
                strncat(hb_at[i2],at1,3);
                if(strlen(hb_seq[i2])>=HBMAX)
                {
                  //  printf("Increase HBMAX\n");
                    printf("ERROR 10\n");
                    exit(10);
                }
            }
        }
    }
    fclose(seqin);
}
/*-------------------------------------------*/
bool solve_eq(float coM[],
              float co1[],
              float co2[],
              float an1,
              float an2,
              float l1,
              float l2,
              float l,
              float res[],
              float res1[],
              int* ret)
/*-------------------------------------------
resolution des eqs:
x x1 + y y1 + z z1 = alpha
x x2 + y y2 + z z2 = beta
x^2 + y^2 + z^2 = l^2
 -------------------------------------------*/
{
    int i;
    float x1,y1,z1,x2,y2,z2,a1,a2,a3,a,b,c,del,rac1,rac2;
    float alpha,beta;

    *ret=true;

    x1=co1[0]-coM[0];
    y1=co1[1]-coM[1];
    z1=co1[2]-coM[2];
    x2=co2[0];
    y2=co2[1];
    z2=co2[2];
    a1=x1*y2-x2*y1;
    a2=y1*z2-y2*z1;
    a3=z1*x2-z2*x1;
    a=a1*a1+a2*a2+a3*a3;
    if(a==0)
    {
      //  printf("a=%f, I cannot do anything\n",a);
        printf("ERROR 11\n");
        exit(11);
    }

    alpha=(float)cos((double)an1)*l1*l;
    beta=(float)cos((double)an2)*l2*l;

    if(a2!=0)
    {
        b=(alpha*z2-beta*z1)*a3-(alpha*y2-beta*y1)*a1;
        c=(alpha*z2-beta*z1)*(alpha*z2-beta*z1)+
            (alpha*y2-beta*y1)*(alpha*y2-beta*y1)-l*l*a2*a2;

        del=b*b-a*c;
        if(del<0)
        {
          //  printf("Error, Delta=%f <0\n",del);
            *ret=false;
            return(false);
        }

        rac1=(-b/a+(float)sqrt((double)del/(a*a)));
        rac2=(-b/a-(float)sqrt((double)del/(a*a)));

        res[0]=rac1;
        res1[0]=rac2;
        res[1]=(alpha*z2-beta*z1+res[0]*a3)/a2;
        res1[1]=(alpha*z2-beta*z1+res1[0]*a3)/a2;
        if(z2!=0.0)
        {
            res[2]=(beta-res[0]*x2-res[1]*y2)/z2;
            res1[2]=(beta-res1[0]*x2-res1[1]*y2)/z2;
        }
        else if(z1!=0.0)
        {
            res[2]=(alpha-res[0]*x1-res[1]*y1)/z1;
            res1[2]=(alpha-res1[0]*x1-res1[1]*y1)/z1;
        }
        else
        {
            printf("ERROR 12\n");
            exit(12);
        }

        for(i=0;i<3;i++)
        {
            res[i]+=coM[i];
            res1[i]+=coM[i];
        }
    }
    else if(a1!=0)
    {
        b=(alpha*y2-beta*y1)*a2-(alpha*x2-beta*x1)*a3;
        c=(alpha*y2-beta*y1)*(alpha*y2-beta*y1)+
            (alpha*x2-beta*x1)*(alpha*x2-beta*x1)-l*l*a1*a1;

        del=b*b-a*c;
        if(del<0)
        {
          //  printf("Error, Delta=%f <0\n",del);
            *ret=false;
            return(false);
        }

        rac1=(-b/a+(float)sqrt((double)del/(a*a)));
        rac2=(-b/a-(float)sqrt((double)del/(a*a)));

        res[3]=rac1;
        res1[3]=rac2;
        res[0]=(alpha*y2-beta*y1+res[2]*a2)/a1;
        res1[0]=(alpha*y2-beta*y1+res1[2]*a2)/a1;
        if(y2!=0.0)
        {
            res[1]=(beta-res[0]*x2-res[2]*z2)/y2;
            res1[1]=(beta-res1[0]*x2-res1[2]*z2)/y2;
        }
        else if(y1!=0.0)
        {
            res[1]=(alpha-res[0]*x1-res[2]*z1)/y1;
            res1[1]=(alpha-res1[0]*x1-res1[2]*z1)/y1;
        }
        else
        {
            printf("ERROR 13\n");
            exit(13);
        }

        for(i=0;i<3;i++)
        {
            res[i]+=coM[i];
            res1[i]+=coM[i];
        }
    }
    else if(a3!=0)
    {
        b=(alpha*x2-beta*x1)*a1-(alpha*z2-beta*z1)*a2;
        c=(alpha*z2-beta*z1)*(alpha*z2-beta*z1)+
            (alpha*x2-beta*x1)*(alpha*x2-beta*x1)-l*l*a3*a3;

        del=b*b-a*c;
        if(del<0)
        {
          //  printf("Error, Delta=%f <0\n",del);
            *ret=false;
            return(false);
        }

        rac1=(-b/a+(float)sqrt((double)del/(a*a)));
        rac2=(-b/a-(float)sqrt((double)del/(a*a)));

        res[1]=rac1;
        res1[1]=rac2;
        res[2]=(alpha*x2-beta*x1+res[1]*a1)/a3;
        res1[2]=(alpha*x2-beta*x1+res1[1]*a1)/a3;
        if(x2!=0.0)
        {
            res[0]=(beta-res[1]*y2-res[2]*z2)/x2;
            res1[0]=(beta-res1[1]*y2-res1[2]*z2)/x2;
        }
        else if(x1!=0.0)
        {
            res[0]=(alpha-res[1]*y1-res[2]*z1)/x1;
            res1[0]=(alpha-res1[1]*y1-res1[2]*z1)/x1;
        }
        else
        {
            printf("ERROR 14\n");
            exit(14);
        }

        for(i=0;i<3;i++)
        {
            res[i]+=coM[i];
            res1[i]+=coM[i];
        }
    }
    return(true);
}
/*-------------------------------------------*/
void allocate_res(int nres,
                  char **seq,
                  char **longseq,
                  char **flagADN,
                  char **hist,
                  char **prefixes,
                  char (**hb_seq)[HBMAX],
                  char (**hb_at)[7*HBMAX],
                  char (**hb_hist)[6*HBMAX],
                  float (**coord)[2][3],
                  float (**coordN)[2][3],
                  float (**coordtot)[9][3])
/*-------------------------------------------*/
{   int i,j,k;
    
    if((*seq =  malloc((nres+1) * sizeof(char)))==NULL)
    {
        printf("Memory not allocated for seq.\n");
        exit(0);
    }
    memset(*seq, '\0', (nres+1) * sizeof(char));
    
    if((*longseq=  malloc((3*nres+1) * sizeof(char)))==NULL)
    {
        printf("Memory not allocated for seq.\n");
        exit(0);
    }
    memset(*longseq, '\0', (3*nres+1) * sizeof(char));
    
    if((*flagADN =  malloc((nres+1) * sizeof(char)))==NULL)
    {
        printf("Memory not allocated for flagADN.\n");
        exit(0);
    }
    memset(*flagADN, '3', nres * sizeof(char));
    
    (*flagADN)[nres]='\0';
    
    if((*hist = malloc((6*nres+1) * sizeof(char)))==NULL)
    {
        printf("Memory not allocated for hist.\n");
        exit(0);
    }
    memset(*hist, '\0', (6*nres+1) * sizeof(char));

    if((*prefixes = malloc((nres+1) * sizeof(char)))==NULL)
    {
        printf("Memory not allocated for prefixes.\n");
        exit(0);
    }
    memset(*prefixes, '\0', (nres+1) * sizeof(char));
        
    if((*hb_seq =  malloc((nres+1) * sizeof(**hb_seq)))==NULL)
    {
        printf("Memory not allocated for hb_seq.\n");
        exit(0);
    }
    memset(*hb_seq, '\0', (nres+1) * sizeof(**hb_seq));
    
    if((*hb_at =  malloc((nres+1) * sizeof(**hb_at)))==NULL)
    {
        printf("Memory not allocated for hb_at.\n");
        exit(0);
    }
    memset(*hb_at, '\0', (nres+1) * sizeof(**hb_at));
    
    if((*hb_hist =  malloc((nres+1) * sizeof(**hb_hist)))==NULL)
    {
        printf("Memory not allocated for hb_hist.\n");
        exit(0);
    }
    memset(*hb_hist, '\0', (nres+1) * sizeof(**hb_hist));
    
    if((*coord = malloc(nres * sizeof(**coord)))==NULL)
    {
        printf("Memory not allocated for coord.\n");
        exit(0);
    }
    for(i=0;i<nres;i++)
        for(j=0;j<2;j++)
            for(k=0;k<3;k++)
                (*coord)[i][j][k]=NONEXIST;
    
    if((*coordN =  malloc(nres * sizeof(**coordN)))==NULL)
    {
        printf("Memory not allocated for coordN.\n");
        exit(0);
    }
    for(i=0;i<nres;i++)
        for(j=0;j<2;j++)
            for(k=0;k<3;k++)
                (*coordN)[i][j][k]=NONEXIST;
    
    if((*coordtot = malloc(nres * sizeof(**coordtot)))==NULL)
    {
        printf("Memory not allocated for coordtot.\n");
        exit(0);
    }
    for(i=0;i<nres;i++)
        for(j=0;j<9;j++)
            for(k=0;k<3;k++)
                (*coordtot)[i][j][k]=NONEXIST;
}
/*-------------------------------------------*/
void reallocate_res(int ires,
                    char **seq,
                    char **longseq,
                    char **flagADN,
                    char **hist,
                    char **prefixes,
                    char (**hb_seq)[HBMAX],
                    char (**hb_at)[7*HBMAX],
                    char (**hb_hist)[6*HBMAX],
                    float (**coord)[2][3],
                    float (**coordN)[2][3],
                    float (**coordtot)[9][3])
/*-------------------------------------------*/
{   int i,j,k;
    
    if((*seq = realloc(*seq,(ires+1) * sizeof(char)))==NULL)
    {
        printf("Memory not reallocated for seq.\n");
        exit(0);
    }
    memset(*seq, '\0', (ires+1) * sizeof(char));
    
    if((*longseq = realloc(*longseq,(3*ires+1) * sizeof(char)))==NULL)
    {
        printf("Memory not reallocated for longseq.\n");
        exit(0);
    }
    memset(*longseq, '\0', (3*ires+1) * sizeof(char));
    
    if((*flagADN = realloc(*flagADN,(ires+1) * sizeof(char)))==NULL)
    {
        printf("Memory not reallocated for flagADN.\n");
        exit(0);
    }
    memset(*flagADN, '3', ires * sizeof(char));
    (*flagADN)[ires]='\0';
    
    if((*hist = realloc(*hist,(6*ires+1) * sizeof(char)))==NULL)
    {
        printf("Memory not reallocated for hist.\n");
        exit(0);
    }
    memset(*hist, '\0', (6*ires+1) * sizeof(char));

    if((*prefixes =  realloc(*prefixes,(ires+1) * sizeof(char)))==NULL)
    {
        printf("Memory not reallocated for prefixes.\n");
        exit(0);
    }
    memset(*prefixes, '\0', (ires+1) * sizeof(char));
        
    if((*hb_seq =  realloc(*hb_seq,(ires+1) * sizeof(**hb_seq)))==NULL)
    {
        printf("Memory not reallocated for hb_seq.\n");
        exit(0);
    }
    memset(*hb_seq, '\0', (ires+1) * sizeof(**hb_seq));
    
    if((*hb_at =  realloc(*hb_at,(ires+1) * sizeof(**hb_at)))==NULL)
    {
        printf("Memory not reallocated for hb_seq.\n");
        exit(0);
    }
    memset(*hb_at, '\0', (ires+1) * sizeof(**hb_at));
    
    if((*hb_hist =  realloc(*hb_hist,(ires+1) * sizeof(**hb_hist)))==NULL)
    {
        printf("Memory not reallocated for hb_hist.\n");
        exit(0);
    }
    memset(*hb_hist, '\0', (ires+1) * sizeof(**hb_hist));
    
    if((*coord = realloc(*coord,ires * sizeof(**coord)))==NULL)
    {
        printf("Memory not reallocated for coord.\n");
        exit(0);
    }
    for(i=0;i<ires;i++)
        for(j=0;j<2;j++)
            for(k=0;k<3;k++)
                (*coord)[i][j][k]=NONEXIST;
    
    if((*coordN = realloc(*coordN,ires * sizeof(**coordN)))==NULL)
    {
        printf("Memory not reallocated for coord.\n");
        exit(0);
    }
    for(i=0;i<ires;i++)
        for(j=0;j<2;j++)
            for(k=0;k<3;k++)
                (*coordN)[i][j][k]=NONEXIST;
    
    if((*coordtot = realloc(*coordtot,ires * sizeof(**coordtot)))==NULL)
    {
        printf("Memory not reallocated for coord.\n");
        exit(0);
    }
    for(i=0;i<ires;i++)
        for(j=0;j<9;j++)
            for(k=0;k<3;k++)
                (*coordtot)[i][j][k]=NONEXIST;
    
}
/*-------------------------------------------*/
void allocate_mem(int memcount,
                char **memseq1,
                char **memseq2,
                char **memlongseq2,
                char **memflagADN2,
                char **memprefix2,
                char (**memhist1)[7],
                char (**memhist2)[7],
                char (**ligne)[1024],
                char (**lignecsv)[1024],
                char (**memhb_hist)[6*HBMAX],
                char (**memhb_seq)[HBMAX],
                char (**memhb_at)[7*HBMAX],
                int **fl,
                int **fldeja,
                int **repeat,
                float **mema1,
                float **memdis,
                float **memangle)
/*-------------------------------------------*/
{   int i;
    
    if((*memseq1 =  malloc(memcount * sizeof(char)))==NULL)
    {printf("Memory not allocated for memseq1.\n");
        exit(0);
    }
    memset(*memseq1, '\0', memcount * sizeof(char));
    
    if((*memseq2 =  malloc(memcount * sizeof(char)))==NULL)
    {printf("Memory not allocated for memseq2.\n");
        exit(0);
    }
    memset(*memseq2, '\0', memcount * sizeof(char));
    
    if((*memlongseq2 =  malloc(memcount * 3 * sizeof(char)))==NULL)
    {printf("Memory not allocated for memlongseq2.\n");
        exit(0);
    }
    memset(*memlongseq2, '\0', memcount * 3 * sizeof(char));
    
    if((*memflagADN2 =  malloc(memcount * sizeof(char)))==NULL)
    {printf("Memory not allocated for memflagADN2.\n");
        exit(0);
    }
    memset(*memflagADN2, '\0', memcount * sizeof(char));
    
    if((*memprefix2 =  malloc(memcount * sizeof(char)))==NULL)
    {printf("Memory not allocated for memprefix2.\n");
        exit(0);
    }
    memset(*memprefix2, '\0', memcount * sizeof(char));
    
    if((*memhist1 =  malloc(memcount * sizeof(*memhist1)))==NULL)
    {printf("Memory not allocated for memhist1.\n");
        exit(0);
    }
    memset(*memhist1, '\0', memcount * sizeof(*memhist1));
    
    if((*memhist2 =  malloc(memcount * sizeof(**memhist2)))==NULL)
    {printf("Memory not allocated for memhist2.\n");
        exit(0);
    }
    memset(*memhist2, '\0', memcount * sizeof(**memhist2));
    
    if((*ligne =  malloc(memcount * memcount * sizeof(**ligne)))==NULL)
    {
        printf("Memory not allocated for ligne.\n");
        exit(0);
    }
    memset(*ligne, '\0', memcount * memcount * sizeof(**ligne));
    
    if((*lignecsv =  malloc(memcount * memcount * sizeof(**lignecsv)))==NULL)
    {
        printf("Memory not allocated for lignecsv.\n");
        exit(0);
    }
    memset(*lignecsv, '\0', memcount * memcount * sizeof(**lignecsv));

    if((*memhb_hist =  malloc(memcount * sizeof(**memhb_hist)))==NULL)
    {printf("Memory not allocated for memhb_hist.\n");
        exit(0);
    }
    memset(*memhb_hist, '\0', memcount * sizeof(**memhb_hist));
    
    if((*memhb_seq =  malloc(memcount * sizeof(**memhb_seq)))==NULL)
    {printf("Memory not allocated for memhb_seq.\n");
        exit(0);
    }
    memset(*memhb_seq, '\0', memcount * sizeof(**memhb_seq));
    
    if((*memhb_at =  malloc(memcount * sizeof(**memhb_at)))==NULL)
    {printf("Memory not allocated for memhb_at.\n");
        exit(0);
    }
    memset(*memhb_at, '\0', memcount * sizeof(**memhb_at));
    
    if((*fl = malloc(memcount * sizeof(int)))==NULL)
    {printf("Memory not allocated for fl.\n");
        exit(0);
    }
    memset(*fl, 0, memcount * sizeof(int));
    
    if((*fldeja = malloc(memcount * memcount * sizeof(int)))==NULL)
    {printf("Memory not allocated for fldeja.\n");
        exit(0);
    }
    memset(*fldeja, 0, memcount * memcount * sizeof(int));
    
    if((*repeat = malloc(memcount * sizeof(int)))==NULL)
    {printf("Memory not allocated for repeat.\n");
        exit(0);
    }
    memset(*repeat, 0, memcount * sizeof(int));
    
    if((*mema1 = malloc(memcount * sizeof(float)))==NULL)
    {printf("Memory not allocated for mema1.\n");
        exit(0);
    }
    for(i=0;i<memcount;i++)(*mema1)[i]=(+NONEXIST);

    if((*memdis = malloc(memcount * sizeof(float)))==NULL)
    {printf("Memory not allocated for memdis.\n");
        exit(0);
    }
    for(i=0;i<memcount;i++)(*memdis)[i]=(-NONEXIST);

    if((*memangle = malloc(memcount * sizeof(float)))==NULL)
    {printf("Memory not allocated for memangle.\n");
        exit(0);
    }
    for(i=0;i<memcount;i++)(*memangle)[i]=(+NONEXIST);
}
/*-------------------------------------------*/
void reallocate_mem(int memcount,
                int memcountold,
                char **memseq1,
                char **memseq2,
                char **memlongseq2,
                char **memflagADN2,
                char **memprefix2,
                char (**memhist1)[7],
                char (**memhist2)[7],
                char (**ligne)[1024],
                char (**lignecsv)[1024],
                char (**memhb_hist)[6*HBMAX],
                char (**memhb_seq)[HBMAX],
                char (**memhb_at)[7*HBMAX],
                int **fl,
                int **fldeja,
                int **repeat,
                float **mema1,
                float **memdis,
                float **memangle)
/*-------------------------------------------*/
{
    int i,ii,ii1;
    ii=memcount-memcountold;
    ii1=memcount*memcount-memcountold*memcountold;
    
    if((*memseq1 = realloc(*memseq1,memcount * sizeof(char)))==NULL)
    {
        printf("Memory not reallocated for memseq1.\n");
        exit(0);
    }
    memset(&(*memseq1)[memcountold], '\0', ii * sizeof(char));
    
    if((*memseq2 = realloc(*memseq2,memcount * sizeof(char)))==NULL)
    {
        printf("Memory not reallocated for memseq2.\n");
        exit(0);
    }
    memset(&(*memseq2)[memcountold], '\0', ii * sizeof(char));
    
    if((*memlongseq2 = realloc(*memlongseq2,memcount * 3 * sizeof(char)))==NULL)
    {
        printf("Memory not reallocated for memseq2.\n");
        exit(0);
    }
    memset(&(*memlongseq2)[memcountold*3], '\0', ii * 3 * sizeof(char));
    
    if((*memflagADN2 = realloc(*memflagADN2,memcount * sizeof(char)))==NULL)
    {
        printf("Memory not reallocated for memflagADN2.\n");
        exit(0);
    }
    memset(&(*memflagADN2)[memcountold], '\0', ii * sizeof(char));
    
    if((*memprefix2 = realloc(*memprefix2,memcount * sizeof(char)))==NULL)
    {
        printf("Memory not reallocated for memprefix2.\n");
        exit(0);
    }
    memset(&(*memprefix2)[memcountold], '\0', ii * sizeof(char));
    
    if((*memhist1 =  realloc(*memhist1,memcount * sizeof(**memhist1)))==NULL)
    {
        printf("Memory not reallocated for memhist1.\n");
        exit(0);
    }
    memset(&(*memhist1)[memcountold], '\0', ii * sizeof(**memhist1));
    
    if((*memhist2 =  realloc(*memhist2,memcount * sizeof(**memhist2)))==NULL)
    {
        printf("Memory not reallocated for memhist2.\n");
        exit(0);
    }
    memset(&(*memhist2)[memcountold], '\0', ii * sizeof(**memhist2));
    
    if((*ligne =  realloc(*ligne,memcount * memcount * sizeof(**ligne)))==NULL)
    {
        printf("Memory not reallocated for ligne.\n");
        exit(0);
    }
    memset(&(*ligne)[memcountold*memcountold], '\0', ii1 * sizeof(**ligne));
    
    if((*lignecsv =  realloc(*lignecsv,memcount * memcount * sizeof(**lignecsv)))==NULL)
    {
        printf("Memory not reallocated for lignecsv.\n");
        exit(0);
    }
    memset(&(*lignecsv)[memcountold*memcountold], '\0', ii1 * sizeof(**lignecsv));
        
    if((*memhb_hist =  realloc(*memhb_hist,memcount * sizeof(**memhb_hist)))==NULL)
    {
        printf("Memory not reallocated for memhb_hist.\n");
        exit(0);
    }
    memset(&(*memhb_hist)[memcountold], '\0', ii * sizeof(**memhb_hist));
    
    if((*memhb_seq =  realloc(*memhb_seq,memcount * sizeof(**memhb_seq)))==NULL)
    {
        printf("Memory not reallocated for memhb_seq.\n");
        exit(0);
    }
    memset(&(*memhb_seq)[memcountold], '\0', ii * sizeof(**memhb_seq));
    
    if((*memhb_at =  realloc(*memhb_at,memcount * sizeof(**memhb_at)))==NULL)
    {
        printf("Memory not reallocated for memhb_at.\n");
        exit(0);
    }
    memset(&(*memhb_at)[memcountold], '\0', ii * sizeof(**memhb_at));
    
    if((*fl = realloc(*fl,memcount * sizeof(int)))==NULL)
    {
        printf("Memory not reallocated for fl.\n");
        exit(0);
    }
    memset(&(*fl)[memcountold], 0, ii * sizeof(int));
    
    if((*fldeja = realloc(*fldeja,memcount * memcount * sizeof(int)))==NULL)
    {
        printf("Memory not reallocated for fldeja.\n");
        exit(0);
    }
    memset(&(*fldeja)[memcountold*memcountold], 0, ii1 * sizeof(int));
    
    if((*repeat = realloc(*repeat,memcount * sizeof(int)))==NULL)
    {
        printf("Memory not reallocated for repeat.\n");
        exit(0);
    }
    memset(&(*repeat)[memcountold], 0, ii * sizeof(int));
    
    if((*mema1 = realloc(*mema1,memcount * sizeof(float)))==NULL)
    {
        printf("Memory not reallocated for mema1.\n");
        exit(0);
    }
    for(i=0;i<ii;i++)(*mema1)[memcountold+i]=(+NONEXIST);
    
    if((*memdis = realloc(*memdis,memcount * sizeof(float)))==NULL)
    {
        printf("Memory not reallocated for memdis.\n");
        exit(0);
    }
    for(i=0;i<ii;i++)(*memdis)[memcountold+i]=(-NONEXIST);

    if((*memangle = realloc(*memangle,memcount * sizeof(float)))==NULL)
    {
        printf("Memory not reallocated for memangle.\n");
        exit(0);
    }
    for(i=0;i<ii;i++)(*memangle)[memcountold+i]=(+NONEXIST);
}
/*-------------------------------------------*/
int main(void)
/*-------------------------------------------*/
{
    int i,cc,nres,ires,j,k,t;
    int m,mm,mem,jj,kk;
    int icc,vv,memchargaro,memarosulfur;
    int aro,doub;
    int comp,i1,i2;
    int deja;
    int *fldeja,*fl,*repeat,memcount,memcountold;
    int num[5][5],fou,q,rep;
    int iimax,iimax1,iim,iima,iip,iiq,iich,iiar;
    int drap,redundancy=(int)false;
    char repert[1024],string[1024];
    char flhbplus[2],fichin[1024];
    char (*hb_seq)[HBMAX],(*hb_at)[7*HBMAX],(*hb_hist)[6*HBMAX];
    char (*memhb_hist)[6*HBMAX],(*memhb_seq)[HBMAX];
    char (*memhb_at)[7*HBMAX];
    char fichlist[1024],fichin1[1024];
    char fichout[1024],fichoutcsv[1024];
    char fichoutcsv_chst[1024];
    char mot[20],mot1[20],mot11[20];
    char *hist,*seq,*flagADN;
    char *longseq,*prefixes;
    char aux1[7],aux2[7],aux3[7],auxh[7],aux4[8];
    char auxlong[4],auxpref[2],str[1024];
    char *memseq1,*memseq2,*memlongseq2;
    char *memflagADN2,*memprefix2;
    char (*memhist1)[7],(*memhist2)[7];
    char (*ligne)[1024],(*lignecsv)[1024];
    char name[2560],prot[40];
    char path2code[1024],*dircode;
    float d,dm,d1,d2,tt[3],ttN[3];
    float dc5,dc6,ff;
    float *mema1,*memdis,*memangle;
    float amax;
    float poid[30];
    float rmsch[5][5];
    float Amax[2],Amaxpi[2],a,a1;
    float CATradcyl,PiPiradcyl,SPiradcyl;
    float (*coord)[2][3],(*coordN)[2][3],(*coordtot)[9][3];
    FILE *seqin,*seqlist,*seqout,*seqoutcsv;
    FILE *seqoutcsv_chst;
    
    dircode = getenv("DIRPI");
    sprintf(path2code,"%s/code",dircode);

    sprintf(fichin,"%s/parameters",dircode);
    seqin=fopen(fichin,"r");
    if(seqin==NULL)
    {   
        printf("\nFile '%s' not existing\n",fichin);
        printf("I will use the default parameters\n");
    }
    else
    {
        while(fscanf(seqin,"%s %f",mot,&ff)!=EOF)
        {
            if(!strcmp(mot,"CATDmax"))CATDmax=ff;
            else if(!strcmp(mot,"CATangmax1"))CATangmax1=ff;
            else if(!strcmp(mot,"PiPiDmax"))PiPiDmax=ff;
            else if(!strcmp(mot,"PiPiangmax1"))PiPiangmax1=ff;
            else if(!strcmp(mot,"SPiDmax"))SPiDmax=ff;
            else if(!strcmp(mot,"SPiangmax1"))SPiangmax1=ff;
            else if(!strcmp(mot,"redundancy"))redundancy=(int)ff;
        }
        fclose(seqin);
        if(redundancy!=0&&redundancy!=1)
        {
            printf("*** Redundancy is either 1 (TRUE) or 0 (FALSE), not %d ***\n",redundancy);
            printf("ERROR 333\n");
            exit(333);
        }
    }

    printf("\nParameters used:\n");
    printf(" Cation/amino/His-pi: interresidue distance < %3.1f \u212B and cylinder of radius %3.1f x radius of aromatic cycle\n",
             CATDmax,CATangmax1);
    printf(" Pi-pi: interresidue distance < %3.1f \u212B and cylinder of radius %3.1f x radius of aromatic cycle\n",
	     PiPiDmax,PiPiangmax1);
    printf(" Sulfur-pi: interresidue distance < %3.1f \u212B and cylinder of radius %3.1f x radius of aromatic cycle\n",
	     SPiDmax,SPiangmax1);
    if(redundancy==0)
        printf(" No redundancy: do not print pi-interactions already shown in other identical chains\n\n");
    else
        printf(" With redundancy: print all pi-interactions even those already shown in identical chains\n\n");

    CATradcyl=CATangmax1;
    PiPiradcyl=PiPiangmax1;
    SPiradcyl=SPiangmax1;
    
    iim=0;
    iima=0;
    a1=0.0;
    
    for(i=0;i<9*3;i++)poid[i]=1.0;
    
    for(i=0;i<5;i++)
        for(j=0;j<5;j++)
        {
            num[i][j]=0;
            rmsch[i][j]=0.0;
        }
    aux1[6]='\0';
    aux2[6]='\0';
    aux3[6]='\0';
    auxh[6]='\0';
    aux4[7]='\0';
    auxlong[3]='\0';
    auxpref[1]='\0';
    
    strcpy(residus,"RKQNHFYWHCMacgtuxy");
    strcpy(rescharg,"RKQNH");
    strcpy(resaro,"FYWH");
    strcpy(resulfur,"CM");
    strcpy(resbase,"acgtuxy");
    strcpy(ttamines,"ACDEFGHIKLMNPQRSTVWY");
    namines=(int)strlen(ttamines);
    
    strcpy(atbase[0],"N1");
    strcpy(atbase[1],"N3");
    strcpy(atbase[2],"C2");
    strcpy(atbase[3],"C6");
    strcpy(atbase[4],"C4");
    strcpy(atbase[5],"C5");
    strcpy(atbase[6],"C8");
    strcpy(atbase[7],"N7");
    strcpy(atbase[8],"N9");
    
    i=(int)(strchr(resbase,'a')-resbase);
    atbamax[i]=9;
    i=(int)(strchr(resbase,'g')-resbase);
    atbamax[i]=9;
    i=(int)(strchr(resbase,'y')-resbase);
    atbamax[i]=9;
    i=(int)(strchr(resbase,'c')-resbase);
    atbamax[i]=6;
    i=(int)(strchr(resbase,'t')-resbase);
    atbamax[i]=6;
    i=(int)(strchr(resbase,'u')-resbase);
    atbamax[i]=6;
    i=(int)(strchr(resbase,'x')-resbase);
    atbamax[i]=6;
    
    i=(int)(strchr(rescharg,'R')-rescharg);
    strcpy(atcharg[i][0],"NH1");
    strcpy(atcharg[i][1],"NH2");
    strcpy(atcharg[i][2],"CZ");
    strcpy(atcharg[i][3],"NE");
    atchmax[i]=4;
    i=(int)(strchr(rescharg,'K')-rescharg);
    strcpy(atcharg[i][0],"NZ");
    strcpy(atcharg[i][1],"CE");
    strcpy(atcharg[i][2],"H1");
    atchmax[i]=3;
    i=(int)(strchr(rescharg,'Q')-rescharg);
    strcpy(atcharg[i][0],"NE2");
    strcpy(atcharg[i][1],"CD");
    strcpy(atcharg[i][2],"OE1");
    atchmax[i]=3;
    i=(int)(strchr(rescharg,'N')-rescharg);
    strcpy(atcharg[i][0],"ND2");
    strcpy(atcharg[i][1],"CG");
    strcpy(atcharg[i][2],"OD1");
    atchmax[i]=3;
    i=(int)(strchr(rescharg,'H')-rescharg);
    strcpy(atcharg[i][0],"ND1");
    strcpy(atcharg[i][1],"CE1");
    strcpy(atcharg[i][2],"NE2");
    atchmax[i]=3;
    
    i=(int)(strchr(resulfur,'C')-resulfur);
    strcpy(atsulfur[i][0],"SG");
    atsulfurmax[i]=1;
    i=(int)(strchr(resulfur,'M')-resulfur);
    strcpy(atsulfur[i][0],"SD");
    atsulfurmax[i]=1;
    
    i=(int)(strchr(resaro,'H')-resaro);
    strcpy(ataro[i][0],"ND1");
    strcpy(ataro[i][1],"CE1");
    strcpy(ataro[i][2],"NE2");
    strcpy(ataro[i][3],"CD2");
    strcpy(ataro[i][4],"CG");
    atarmax[i]=5;
    i=(int)(strchr(resaro,'Y')-resaro);
    strcpy(ataro[i][0],"CG");
    strcpy(ataro[i][1],"CD1");
    strcpy(ataro[i][2],"CD2");
    strcpy(ataro[i][3],"CE1");
    strcpy(ataro[i][4],"CE2");
    strcpy(ataro[i][5],"CZ");
    atarmax[i]=6;
    i=(int)(strchr(resaro,'F')-resaro);
    strcpy(ataro[i][0],"CG");
    strcpy(ataro[i][1],"CD1");
    strcpy(ataro[i][2],"CD2");
    strcpy(ataro[i][3],"CE1");
    strcpy(ataro[i][4],"CE2");
    strcpy(ataro[i][5],"CZ");
    atarmax[i]=6;
    i=(int)(strchr(resaro,'W')-resaro);
    strcpy(ataro[i][0],"CZ2");
    strcpy(ataro[i][1],"CH2");
    strcpy(ataro[i][2],"CZ3");
    strcpy(ataro[i][3],"CE3");
    strcpy(ataro[i][4],"CD2");
    strcpy(ataro[i][5],"CE2");
    strcpy(ataro[i][6],"NE1");
    strcpy(ataro[i][7],"CD1");
    strcpy(ataro[i][8],"CG");
    atarmax[i]=9;
    
    nres=100;
    allocate_res(nres,&seq,&longseq,&flagADN,&hist,&prefixes,&hb_seq,&hb_at,&hb_hist,&coord,&coordN,&coordtot);
    
    memcount=MEMMAX;
    allocate_mem(memcount,&memseq1,&memseq2,&memlongseq2,&memflagADN2,
                 &memprefix2,&memhist1,&memhist2,&ligne,&lignecsv,&memhb_hist,&memhb_seq,
                 &memhb_at,&fl,&fldeja,&repeat,&mema1,&memdis,&memangle);
    
    printf("Name of directory containing PDB files to be analyzed: ");
    if(scanf("%s",repert)!=1)
    {
        printf("*** Incorrect answer: '%s'\n",repert);
        exit(555);
    }
    i=(int)strlen(repert);
    if(repert[i-1]!='/')
    {
        repert[i]='/';
        repert[i+1]='\0';
    }
    
    dc6=-FLT_MAX;
    dc5=dc6;
    for(i=0;i<2;i++)Amaxpi[i]=dc6;
    
    for(i=0;i<strlen(resbase)-3;i++)
    {
        sprintf(fichin,"%s/optimizedPDBs/%c_opt.pdb",path2code,resbase[i]);
        read_protADN(fichin,&nres,seq,longseq,prefixes,hist,Bcoord[i],BcoordN[i],Bcoordtot[i],flagADN,name);
        if(nres!=1)
        {
            //  printf("Error, more than one nucleobase in %s\n",fichin);
            printf("ERROR 15\n");
            exit(15);
        }
        if(flagADN[0]!='1'&&flagADN[0]!='2')
        {
            //  printf("Error, no DNA in %s\n",fichin);
            printf("ERROR 16\n");
            exit(16);
        }
        
        for(j=0;j<6;j++)
        {
            d=dist2(Bcoordtot[i][0][j],Bcoord[i][0][j]); //A CHECKER !!
            if(d>dc6)dc6=d;
            
            for(k=0;k<3;k++)tt[k]=(Bcoordtot[i][0][j][k]-Bcoord[i][0][0][k])*2.0;
            for(k=0;k<3;k++)ttN[k]=BcoordN[i][0][0][k];
            for(k=0;k<3;k++)ttN[k]*=CATDmax;
            for(k=0;k<3;k++)ttN[k]+=tt[k];
            a=M_PI/2.0-angle(tt,ttN,true);
            if(a>Amaxpi[0])Amaxpi[0]=a;
        }
        
        if(resbase[i]=='a'||resbase[i]=='g')
            for(j=4;j<9;j++)
            {
                d=dist2(Bcoordtot[i][0][j],Bcoord[i][0][1]);
                if(d>dc5)dc5=d;
                
                for(k=0;k<3;k++)tt[k]=(Bcoordtot[i][0][j][k]-Bcoord[i][0][0][k])*2.0;
                for(k=0;k<3;k++)ttN[k]=BcoordN[i][0][0][k];
                for(k=0;k<3;k++)ttN[k]*=CATDmax;
                for(k=0;k<3;k++)ttN[k]+=tt[k];
                a=M_PI/2.0-angle(tt,ttN,true);
                if(a>Amaxpi[1])Amaxpi[1]=a;
            }
    }
    
    for(i=0;i<2;i++)Amax[i]=180.0*Amaxpi[i]/M_PI;
    
    for(i=0;i<strlen(resaro);i++)
    {
        sprintf(fichin,
                "%s/optimizedPDBs/%c_opt.pdb",path2code,resaro[i]);
        read_protADN(fichin,&nres,seq,longseq,prefixes,hist,Acoord[i],AcoordN[i],Acoordtot[i],flagADN,name);
        if(nres!=1)
        {
          //  printf("erreur, + d'1 base dans %s\n",fichin);
            printf("ERROR 17\n");
            exit(17);
        }
        if(flagADN[0]!='0')
        {
          //  printf("Error, no DNA in %s\n",fichin);
            printf("ERROR 18\n");
            exit(18);
        }
    }
    
    for(i=0;i<=strlen(rescharg)+1;i++)
    {
        if(i==strlen(rescharg))
            sprintf(fichin,"%s/optimizedPDBs/H_np_opt_1.pdb",path2code);
        else if(i==strlen(rescharg)+1)
            sprintf(fichin,"%s/optimizedPDBs/H_np_opt_2.pdb",path2code);
        else if(rescharg[i]=='H')
            sprintf(fichin,"%s/optimizedPDBs/H_p_opt.pdb",path2code);
        else
            sprintf(fichin,"%s/optimizedPDBs/%c_opt.pdb",path2code,rescharg[i]);
        read_protADN(fichin,&nres,seq,longseq,prefixes,hist,Ccoord[i],CcoordN[i],Ccoordtot[i],flagADN,name);
        if(nres!=1)
        {
          //  printf("Error, more than one nucleobase in %s\n",fichin);
            printf("ERROR 19\n");
            exit(19);
        }
        if(flagADN[0]!='0')
        {
           // printf("Error, something else than an amino acid in %s\n",fichin);
            printf("ERROR 20\n");
            exit(20);
        }
    }
    
    for(i=0;i<strlen(resulfur);i++)
    {
        sprintf(fichin,
                "%s/optimizedPDBs/%c.pdb",path2code,resulfur[i]);
        read_protADN(fichin,&nres,seq,longseq,prefixes,hist,Dcoord[i],DcoordN[i],Dcoordtot[i],flagADN,name);
        if(nres!=1)
        {
          //  printf("error, more than 1 base in %s\n",fichin);
            printf("ERROR 21\n");
            exit(21);
        }
        if(flagADN[0]!='0')
        {
           // printf("Error, no DNA in %s\n",fichin);
            printf("ERROR 22\n");
            exit(22);
        }
    }
    
    sprintf(fichout,"%s%s",repert,"PInteract.txt");
    seqout=fopen(fichout,"w");
    if(seqout==NULL)
    {
        printf("*** '%s' cannot be opened ***\n",fichout);
        printf("ERROR 23\n");
        exit(23);
    }
        
    sprintf(fichoutcsv,"%s%s",repert,"PInteract.csv");
    seqoutcsv=fopen(fichoutcsv,"w");
    if(seqoutcsv==NULL)
    {
        printf("*** '%s' cannot be opened ***\n",fichoutcsv);
        printf("ERROR 24\n");
        exit(24);
    }
    sprintf(fichoutcsv_chst,"%s%s",repert,"PInteract1.csv");
    seqoutcsv_chst=fopen(fichoutcsv_chst,"w");
    if(seqoutcsv_chst==NULL)
    {
        printf("*** '%s' cannot be opened ***\n",fichoutcsv_chst);
        printf("ERROR 24\n");
        exit(24);
    }
    
    sprintf(string," cd %s ; ls *.hb2*",repert);
    if(system(string)==0)flhbplus[0]=(char)true;
    else flhbplus[0]=(char)false;

    fprintf(seqout,
        "# Cation-pi, His-pi and amino-pi interactions defined by an \n#  interresidue distance D < %3.1f \u212B and a cylinder\n#  with radius %3.1f times the radius of the aromatic cycle.\n# A: angle between the normal to the aromatic cycle and the vector between\n#  the center of the cycle and the net or partial charge (maximum angle)\n# B: angle between the aromatic cycle and the plane of the net or partial\n#  charge (B=0° is a stacked conformation; B=90° is a T-shaped conformation)\n",
        CATDmax,CATradcyl);
    fprintf(seqout,
        "#\n# Pi-pi interactions defined by an \n#  interresidue distance D < %3.1f \u212B and a cylinder\n#  with radius %3.1f times the radius of the first aromatic cycle.\n# A: angle between the normal to the first aromatic cycle and the\n#  vector between the center of the two aromatic cycles (maximum angle)\n# B: angle between the two aromatic cycles\n#  (B=0° is a stacked conformation; B=90° is a T-shaped conformation)\n",
        PiPiDmax,PiPiradcyl);
    fprintf(seqout,
        "#\n# Sulfur-pi interactions defined by an \n#  interresidue distance D < %3.1f \u212B and a cylinder\n#  with radius %3.1f times the radius of the aromatic cycle.\n# A: angle between the normal to the aromatic cycle and the vector between\n#  the center of the cycle and the sulfur atom (maximum angle)\n",
        SPiDmax,SPiradcyl);
    if(flhbplus[0]==(char)true)
        fprintf(seqout,
            "#\n# H-bonds computed by HBplus; HB1 is with respect to residue/nucleobase 1;\n#  HB2 is with respect to residue/nucleobase 2\n");
    if(redundancy==0)
        fprintf(seqout,
            "#\n# No redundancy: do not print pi-interactions already shown in other identical chains\n");
    else
        fprintf(seqout,
            "#\n# With redundancy: print all pi-interactions even those already shown in identical chains\n");
    fprintf(seqout,
        "#---------------------------------------------------------------\n");
        
    fprintf(seqoutcsv,
        "PDB_id,Type,TypeInt,Seq1,Nr1,LigName,Seq2,Cycle,Nr2,Dist(\u212B),At1,At2,A(°),Amax(°),B(°)");
    if(flhbplus[0]==(char)true)
    {
        for(i=1;i<=3;i++)
        {
            fprintf(seqoutcsv,",HB1_%d_Seq,HB1_%d_Nr,HB1_%d_At1,HB1_%d_At",i,i,i,i);
        }
        for(i=1;i<=3;i++)
        {
            fprintf(seqoutcsv,",HB2_%d_Seq,HB2_%d_Nr,HB2_%d_At2,HB2_%d_At",i,i,i,i);
        }
    }
    
    fprintf(seqoutcsv,"\n");
    fprintf(seqoutcsv_chst,"PDB_id,Type,Seq1,Nr1,Seq2,Nr2,Seq3,Nr3,At2,At3\n");

    
    sprintf(string," cd %s ; ls *.pdb* > list",repert);
    if(system(string)!=0)
    {
        printf("*** Command '%s' failed\n",string);
        printf("ERROR 45\n");
        exit(45);
    }
    sprintf(fichlist,"%slist",repert);
    seqlist=fopen(fichlist,"r");
    if(seqlist==NULL)
    {
        printf("*** '%s' missing\n",fichlist);
        printf("ERROR 25\n");
        exit(25);
    }

    cc=0;
    while(fscanf(seqlist,"%s",fichin)!=EOF)
    {
        if(!strcmp(&fichin[strlen(fichin)-3],".gz"))
        {
            sprintf(string," gunzip %s%s",repert,fichin);
            if(system(string)!=0)
            {
                printf("*** Command '%s' failed\n",string);
                printf("ERROR 46\n");
                exit(46);
            }
            fichin[strlen(fichin)-3]='\0';
        }
        if(!strcmp(fichin,"pdb"))jj=3;
        else jj=0;
        i=(int)(strchr(fichin,'.')-fichin);
        if(i>=0&&i<strlen(fichin))kk=i;
        else kk=(int)strlen(fichin);
        strncpy(prot,&fichin[jj],kk-jj);
        prot[kk-jj]='\0';
        printf("%s",prot);
        sprintf(fichin1,"%s%s",repert,fichin);
        read_protADN_length(fichin1,&ires);

        reallocate_res(ires,&seq,&longseq,&flagADN,&hist,&prefixes,&hb_seq,&hb_at,&hb_hist,&coord,&coordN,&coordtot);
        
        read_protADN(fichin1,&nres,seq,longseq,prefixes,hist,coord,coordN,coordtot,flagADN,name);

        if(nres > ires)
        {
            printf("error, ires=%d < nres=%d\n",ires,nres);
            printf("%s\n",seq);
            printf("ERROR 26\n");
            exit(26);
        }
        printf(" - %d residues\n",nres);
        
        vv=(int)strlen(seq);
        if((int)strlen(longseq)!=vv*3)
        {
           // printf("error, strlen(seq)=%d and strlen(longseq)=%d\n",vv,(int)strlen(longseq));
            printf("ERROR 27\n");
            exit(27);
        }
        if(vv!=nres)
        {
          //  printf("Error, inconsistent numbers of residues: %d and %d\n",vv,nres);
            printf("ERROR 28\n");
            exit(28);
        }
        icc=0;
        for(m=0;m<nres;m++)
        {
            if(flagADN[m]!='2')continue;
            icc++;
            if(icc>=200)
            {
                //  printf("Increase 200\n");
                printf("ERROR 29\n");
                exit(29);
            }
        }
        
        read_hbonds(fichin,repert,nres,seq,hist,hb_seq,hb_hist,hb_at,flhbplus);
        
//        sprintf(string,"gzip %s",fichin1);
//        if(system(string)!=0)
//        {
//            printf("*** Command '%s' failed\n",string);
//            printf("ERROR 45\n");
//            exit(45);
//        }

        mem=0;
        for(i=0;i<nres;i++)
        {
            if(flagADN[i]>='1')continue;
            iich=(int)(strchr(rescharg,seq[i])-rescharg);
            if(iich<0||iich>strlen(rescharg))continue;
            iimax=atchmax[iich];
            if(seq[i]=='K'||seq[i]=='Q'||seq[i]=='N')iimax--;
            for(j=0;j<nres;j++)
            {
                if(i==j)continue;
                if((seq[i]=='H'&&seq[j]=='H'))continue;
                iiar=(int)(strchr(resaro,seq[j])-resaro);
                if(iiar>=0&&iiar<strlen(resaro))
                {
                    iimax1=atarmax[iiar];
                    aro=true;
                    if(seq[j]=='W')doub=true;
                    else doub=false;
                }
                else
                {
                    iiar=(int)(strchr(resbase,seq[j])-resbase);
                    if(iiar>=0&&iiar<strlen(resbase))
                    {
                        iimax1=atbamax[iiar];
                        aro=false;
                        if(seq[j]=='a'||seq[j]=='g'||seq[j]=='y')doub=true;
                        else doub=false;
                    }
                    else continue;
                }
                d=(float)NONEXIST;
                for(iip=0;iip<iimax;iip++)
                {
                    for(iiq=0;iiq<iimax1;iiq++)
                    {
                        dm=dist2(coordtot[i][iip],coordtot[j][iiq]);
                        if(dm<d)
                        {
                            d=dm;
                            iim=iip;
                            iima=iiq;
                        }
                    }
                }
                if(!doub)q=0;
                else
                {
                    d1=dist2(coordtot[i][iim],coord[j][0]);
                    d2=dist2(coordtot[i][iim],coord[j][1]);
                    if(d1<d2)q=0;
                    else q=1;
                }
                
                if(d<=(CATDmax*CATDmax))
                {
                    for(k=0;k<3;k++)ttN[k]=coordN[j][q][k];
                    if(seq[i]=='H')
                        for(k=0;k<3;k++)tt[k]=coord[i][0][k]-coord[j][q][k];
                    else
                        for(k=0;k<3;k++)tt[k]=coordtot[i][iim][k]-coord[j][q][k];
                    d1=dist2(tt,ttN);
                    for(k=0;k<3;k++)tt[k]=(-tt[k]);
                    d2=dist2(tt,ttN);
                    if(d2>d1)for(k=0;k<3;k++)tt[k]=(-tt[k]);
                    a=angle(tt,ttN,true);
                    
                    if(aro)anglemaxaro(d,&amax,coordtot[j][iima],coord[j][q],coordN[j][q],CATradcyl);
                    else anglemax(iiar,q,d,&amax,CATradcyl);
                    
                    if(a<=amax)
                    {
                        if(coordN[i][0][0]!=NONEXIST&&coordN[i][0][1]!=NONEXIST&&coordN[i][0][2]!=NONEXIST)
                        {
                            for(k=0;k<3;k++)tt[k]=coordN[i][0][k];
                            d1=dist2(tt,ttN);
                            for(k=0;k<3;k++)tt[k]=(-tt[k]);
                            d2=dist2(tt,ttN);
                            if(d2>d1)for(k=0;k<3;k++)tt[k]=(-tt[k]);
                            a1=angle(tt,ttN,true);
                            a1=a1*180.0/M_PI;
                        }
                        else a1=NONEXIST;
                        
                        strncpy(aux1,&hist[6*i],6);
                        strncpy(aux2,&hist[6*j],6);
                        
                        for(m=0;m<mem;m++)
                        {
                            if(seq[i]==memseq1[m]&&seq[j]==memseq2[m]&&
                               !strncmp(&aux1[1],&memhist1[m][1],5)&&
                               !strncmp(&aux2[1],&memhist2[m][1],5))break;
                        }
                        if(redundancy==0&&m<mem)continue;
                        if(flagADN[j]=='2')
                        {
                            for(m=0;m<mem;m++)
                            {
                                if(memflagADN2[m]=='2'&&seq[i]==memseq1[m]&&seq[j]==memseq2[m]&&
                                   !strncmp(&aux1[1],&memhist1[m][1],5)&&!strncmp(&memlongseq2[m*3],&longseq[j*3],3))break;
                            }
                            if(redundancy==0&&m<mem)continue;
                        }
                        
                        rep=false;
                        d=(float)sqrt((double)d);
                        a=a*180.0/M_PI;
                        amax=amax*180.0/M_PI;
                        fou=1;
                        
                        strcpy(memhist1[mem],aux1);
                        strcpy(memhist2[mem],aux2);
                        memseq1[mem]=seq[i];
                        memseq2[mem]=seq[j];
                        strncpy(&memlongseq2[mem*3],&longseq[j*3],3);
                        memprefix2[mem]=prefixes[j];
                        strncpy(auxlong,&longseq[j*3],3);
                        if(prefixes[j]=='\0')auxpref[0]=' ';
                        else auxpref[0]=prefixes[j];
                        memflagADN2[mem]=flagADN[j];
                        if(seq[i]!='K')mema1[mem]=a1;
                        repeat[mem]=rep;
                        
                        memdis[mem]=d;
                        memangle[mem]=a;
                        fl[mem]=fou;
                        
                        if(seq[j]=='H'||q==1)strcpy(mot,"C5");
                        else strcpy(mot,"C6");
                        if(seq[i]!='K')
                        {
                            sprintf(mot1,"B=%5.1f°",a1);
                            sprintf(mot11,"%.1f",a1);
                        }
                        else
                        {
                            strcpy(mot1,"        ");
                            strcpy(mot11,"");
                        }
                        
                        k=(int)(strchr(rescharg,seq[i])-rescharg);
                        strcpy(aux3,atcharg[k][iim]);
                        if(flagADN[j]=='0')
                        {
                            k=(int)(strchr(resaro,seq[j])-resaro);
                            strcpy(aux4,ataro[k][iima]);
                        }
                        else strcpy(aux4,atbase[iima]);
                        
                        if(flagADN[j]=='2')
                        {
                            sprintf(ligne[mem],"  %c (%s) - %s %c [%s] (%s) : D=%4.2f\u212B {%3s-%3s} A=%4.1f° (%4.1f°) %s ",
                                    seq[i],aux1,auxlong,seq[j],mot,aux2,d,aux3,aux4,a,amax,mot1);
                            if(iich<2)sprintf(lignecsv[mem],"cation-pi,");
                            else if(iich<4)sprintf(lignecsv[mem],"amino-pi,");
                            else sprintf(lignecsv[mem],"His-pi,");
                            sprintf(lignecsv[mem],"%s%c,%s,%s,%c,%s,%s,%.2f,%s,%s,%.1f,%.1f,%s",
                                    lignecsv[mem],seq[i],aux1,auxlong,seq[j],mot,aux2,d,aux3,aux4,a,amax,mot11);
                        }
                        else
                        {
                            sprintf(ligne[mem],"  %c (%s) - %c [%s] (%s) : D=%4.2f\u212B {%3s-%3s} A=%4.1f° (%4.1f°) %s ",
                                    seq[i],aux1,seq[j],mot,aux2,d,aux3,aux4,a,amax,mot1);
                            if(iich<2)sprintf(lignecsv[mem],"cation-pi,");
                            else if(iich<4)sprintf(lignecsv[mem],"amino-pi,");
                            else sprintf(lignecsv[mem],"His-pi,");
                            sprintf(lignecsv[mem],"%s%c,%s,,%c,%s,%s,%.2f,%s,%s,%.1f,%.1f,%s",
                                    lignecsv[mem],seq[i],aux1,seq[j],mot,aux2,d,aux3,aux4,a,amax,mot11);
                        }
                        strcpy(memhb_seq[mem],hb_seq[i]);
                        strcpy(memhb_hist[mem],hb_hist[i]);
                        strcpy(memhb_at[mem],hb_at[i]);
                        for(k=0;k<strlen(hb_seq[i]);k++)
                        {
                            strncpy(aux3,&hb_hist[i][k*6],6);
                            strncpy(aux4,&hb_at[i][k*7],7);
                            if(k==0)sprintf(str," HB1: %c {%s} (%s)",hb_seq[i][k],aux4,aux3);
                            else sprintf(str,"; %c {%s} (%s)",hb_seq[i][k],aux4,aux3);
                            strcat(ligne[mem],str);
                            aux4[3]='\0';
                            sprintf(str,",%c,%s,%s,%3s",hb_seq[i][k],aux3,&aux4[0],&aux4[4]);
                            strcat(lignecsv[mem],str);
                        }
                        if(strlen(hb_seq[i])>3)
                            printf("*** WARNING *** %c (%s) in %s makes more than 3 H-bonds\n",seq[i],aux1,prot);
                        if(flhbplus[0]==(char)true)
                        {
                            for(k=strlen(hb_seq[i]);k<3;k++)
                            {
                                sprintf(str,",,,,");
                                strcat(lignecsv[mem],str);
                            }
                        }
                        for(k=0;k<strlen(hb_seq[j]);k++)
                        {
                            strncpy(aux3,&hb_hist[j][k*6],6);
                            strncpy(aux4,&hb_at[j][k*7],7);
                            if(k==0)sprintf(str," HB2: %c {%s} (%s)",hb_seq[j][k],aux4,aux3);
                            else sprintf(str,"; %c {%s} (%s)",hb_seq[j][k],aux4,aux3);
                            strcat(ligne[mem],str);
                            aux4[3]='\0';
                            sprintf(str,",%c,%s,%s,%s",hb_seq[j][k],aux3,&aux4[0],&aux4[4]);
                            strcat(lignecsv[mem],str);
                        }
                        if(strlen(hb_seq[j])>3)
                            printf("*** WARNING *** %c (%s) in %s makes more than 3 H-bonds\n",seq[j],aux2,prot);
                        if(flhbplus[0]==(char)true)
                        {
                            for(k=strlen(hb_seq[j]);k<3;k++)
                            {
                                sprintf(str,",,,,");
                                strcat(lignecsv[mem],str);
                            }
                        }
                        
                        strcat(ligne[mem],"\n");
                        strcat(lignecsv[mem],"\n");
                        mem++;
                        if(mem>=memcount)
                        {
                            memcountold=memcount;
                            memcount+=MEMMAX;
                            reallocate_mem(memcount,memcountold,&memseq1,&memseq2,&memlongseq2,&memflagADN2,
                                           &memprefix2,&memhist1,&memhist2,&ligne,&lignecsv,&memhb_hist,&memhb_seq,
                                           &memhb_at,&fl,&fldeja,&repeat,&mema1,&memdis,&memangle);
                        }
                    }
                }
            }
        }

        memchargaro=mem;
        for(i=0;i<nres;i++)
        {
            if(flagADN[i]>='1')continue;
            iich=(int)(strchr(resaro,seq[i])-resaro);
            if(iich<0||iich>strlen(resaro))continue;
            iimax=atarmax[iich];
            for(j=0;j<nres;j++)
            {
                iiar=(int)(strchr(resaro,seq[j])-resaro);
                if(iiar>=0&&iiar<strlen(resaro)&&j>i)
                {
                    iimax1=atarmax[iiar];
                    aro=true;
                    if(seq[j]=='W')doub=true;
                    else doub=false;
                }
                else
                {
                    iiar=(int)(strchr(resbase,seq[j])-resbase);
                    if(iiar>=0&&iiar<strlen(resbase))
                    {
                        iimax1=atbamax[iiar];
                        aro=false;
                        if(seq[j]=='a'||seq[j]=='g'||seq[j]=='y')doub=true;
                        else doub=false;
                    }
                    else continue;
                }
                d=(float)NONEXIST;
                for(iip=0;iip<iimax;iip++)
                {
                    for(iiq=0;iiq<iimax1;iiq++)
                    {
                        dm=dist2(coordtot[i][iip],coordtot[j][iiq]);
                        if(dm<d)
                        {
                            d=dm;
                            iim=iip;
                            iima=iiq;
                        }
                    }
                }
                if(!doub)q=0;
                else
                {
                    d1=dist2(coordtot[i][iim],coord[j][0]);
                    d2=dist2(coordtot[i][iim],coord[j][1]);
                    if(d1<d2)q=0;
                    else q=1;
                }
                if(d<=(PiPiDmax*PiPiDmax))
                {
                    for(k=0;k<3;k++)ttN[k]=coordN[j][q][k];
                    for(k=0;k<3;k++)tt[k]=coordtot[i][iim][k]-coordtot[j][iima][k];
                    d1=dist2(tt,ttN);
                    for(k=0;k<3;k++)tt[k]=(-tt[k]);
                    d2=dist2(tt,ttN);
                    if(d2>d1)for(k=0;k<3;k++)tt[k]=(-tt[k]);
                    a=angle(tt,ttN,true);
                    
                    if(aro)
                        anglemaxaro(d,&amax,coordtot[j][iima],coord[j][q],coordN[j][q],PiPiradcyl);
                    else anglemax(iiar,q,d,&amax,PiPiradcyl);
                    
                    if(a<=amax)
                    {
                        if(coordN[i][0][0]!=NONEXIST&&coordN[i][0][1]!=NONEXIST&&
                           coordN[i][0][2]!=NONEXIST)
                        {
                            for(k=0;k<3;k++)tt[k]=coordN[i][0][k];
                            d1=dist2(tt,ttN);
                            for(k=0;k<3;k++)tt[k]=(-tt[k]);
                            d2=dist2(tt,ttN);
                            if(d2>d1)for(k=0;k<3;k++)tt[k]=(-tt[k]);
                            a1=angle(tt,ttN,true);
                            a1=a1*180.0/M_PI;
                        }
                        else a1=NONEXIST;
                        
                        strncpy(aux1,&hist[6*i],6);
                        strncpy(aux2,&hist[6*j],6);
                        
                        for(m=memchargaro;m<mem;m++)
                        {
                            if(seq[i]==memseq1[m]&&seq[j]==memseq2[m]&&
                               !strncmp(&aux1[1],&memhist1[m][1],5)&&
                               !strncmp(&aux2[1],&memhist2[m][1],5))break;
                            if(seq[i]==memseq2[m]&&seq[j]==memseq1[m]&&
                               !strncmp(&aux1[1],&memhist2[m][1],5)&&
                               !strncmp(&aux2[1],&memhist1[m][1],5))break;
                        }
                        if(redundancy==0&&m<mem)continue;
                        if(flagADN[j]=='2')
                        {
                            for(m=memchargaro;m<mem;m++)
                            {
                                if(memflagADN2[m]=='2'&&
                                   seq[i]==memseq1[m]&&seq[j]==memseq2[m]&&
                                   !strncmp(&aux1[1],&memhist1[m][1],5)&&
                                   !strncmp(&memlongseq2[m*3],&longseq[j*3],3))break;
                            }
                            if(redundancy==0&&m<mem)continue;
                        }
                        rep=false;
                        d=(float)sqrt((double)d);
                        a=a*180.0/M_PI;
                        amax=amax*180.0/M_PI;
                        fou=1;
                        
                        strcpy(memhist1[mem],aux1);
                        strcpy(memhist2[mem],aux2);
                        memseq1[mem]=seq[i];
                        memseq2[mem]=seq[j];
                        strncpy(&memlongseq2[mem*3],&longseq[j*3],3);
                        memprefix2[mem]=prefixes[j];
                        strncpy(auxlong,&longseq[j*3],3);
                        if(prefixes[j]=='\0')auxpref[0]=' ';
                        else auxpref[0]=prefixes[j];
                        memflagADN2[mem]=flagADN[j];
                        if(seq[i]!='K'&&seq[j]!='K')mema1[mem]=a1;
                        repeat[mem]=rep;
                        
                        memdis[mem]=d;
                        memangle[mem]=a;
                        fl[mem]=fou;
                        
                        if(seq[j]=='H'||q==1)strcpy(mot,"C5");
                        else strcpy(mot,"C6");
                        sprintf(mot1,"B=%5.1f°",a1);
                        sprintf(mot11,"%.1f",a1);
                        
                        k=(int)(strchr(resaro,seq[i])-resaro);
                        strcpy(aux3,ataro[k][iim]);
                        if(flagADN[j]=='0')
                        {
                            k=(int)(strchr(resaro,seq[j])-resaro);
                            strcpy(aux4,ataro[k][iima]);
                        }
                        else strcpy(aux4,atbase[iima]);
                        
                        if(flagADN[j]=='2')
                        {
                            sprintf(ligne[mem],
                                    "  %c (%s) - %s %c (%s) : D=%4.2f\u212B {%3s-%3s} A=%4.1f° (%4.1f°) %s ",
                                    seq[i],aux1,auxlong,seq[j],aux2,d,aux3,aux4,a,amax,mot1);
                            sprintf(lignecsv[mem],"%s,%c,%s,%s,%c,,%s,%.2f,%s,%s,%.1f,%.1f,%s",
                                    "pi-pi",seq[i],aux1,auxlong,seq[j],aux2,d,aux3,aux4,a,amax,mot11);
                        }
                        else
                        {
                            sprintf(ligne[mem],
                                    "  %c (%s) - %c (%s) : D=%4.2f\u212B {%3s-%3s} A=%4.1f° (%4.1f°) %s ",
                                    seq[i],aux1,seq[j],aux2,d,aux3,aux4,a,amax,mot1);
                            sprintf(lignecsv[mem],"%s,%c,%s,,%c,,%s,%.2f,%s,%s,%.1f,%.1f,%s",
                                    "pi-pi",seq[i],aux1,seq[j],aux2,d,aux3,aux4,a,amax,mot11);
                        }
                        strcpy(memhb_seq[mem],hb_seq[i]);
                        strcpy(memhb_hist[mem],hb_hist[i]);
                        strcpy(memhb_at[mem],hb_at[i]);
                        for(k=0;k<strlen(hb_seq[i]);k++)
                        {
                            strncpy(aux3,&hb_hist[i][k*6],6);
                            strncpy(aux4,&hb_at[i][k*7],7);
                            if(k==0)sprintf(str," HB1: %c {%s} (%s)",hb_seq[i][k],aux4,aux3);
                            else sprintf(str,"; %c {%s} (%s)",hb_seq[i][k],aux4,aux3);
                            strcat(ligne[mem],str);
                            
                            aux4[3]='\0';
                            sprintf(str,",%c,%s,%s,%s",hb_seq[i][k],aux3,&aux4[0],&aux4[4]);
                            strcat(lignecsv[mem],str);
                        }
                        if(strlen(hb_seq[i])>3)
                            printf("*** WARNING *** %c (%s) in %s makes more than 3 H-bonds\n",seq[i],aux1,prot);
                        if(flhbplus[0]==(char)true)
                        {
                            for(k=strlen(hb_seq[i]);k<3;k++)
                            {
                                sprintf(str,",,,,");
                                strcat(lignecsv[mem],str);
                            }
                        }
                        for(k=0;k<strlen(hb_seq[j]);k++)
                        {
                            strncpy(aux3,&hb_hist[j][k*6],6);
                            strncpy(aux4,&hb_at[j][k*7],7);
                            if(k==0)sprintf(str," HB2: %c {%s} (%s)",hb_seq[j][k],aux4,aux3);
                            else sprintf(str,"; %c {%s} (%s)",hb_seq[j][k],aux4,aux3);
                            strcat(ligne[mem],str);
                            aux4[3]='\0';
                            sprintf(str,",%c,%s,%s,%s",hb_seq[j][k],aux3,&aux4[0],&aux4[4]);
                            strcat(lignecsv[mem],str);
                        }
                        if(strlen(hb_seq[j])>3)
                            printf("*** WARNING *** %c (%s) in %s makes more than 3 H-bonds\n",seq[j],aux2,prot);
                        if(flhbplus[0]==(char)true)
                        {
                            for(k=strlen(hb_seq[j]);k<3;k++)
                            {
                                sprintf(str,",,,,");
                                strcat(lignecsv[mem],str);
                            }
                        }
                        
                        strcat(ligne[mem],"\n");
                        strcat(lignecsv[mem],"\n");
                        mem++;
                        if(mem>=memcount)
                        {
                            memcountold=memcount;
                            memcount+=MEMMAX;
                            reallocate_mem(memcount,memcountold,&memseq1,&memseq2,&memlongseq2,&memflagADN2,
                                           &memprefix2,&memhist1,&memhist2,&ligne,&lignecsv,&memhb_hist,
                                           &memhb_seq,&memhb_at,&fl,&fldeja,&repeat,&mema1,&memdis,&memangle);
                        }
                    }
                }
            }
        }
        

        memarosulfur=mem;
        for(i=0;i<nres;i++)
        {
            if(flagADN[i]>='1')continue;
            iich=(int)(strchr(resulfur,seq[i])-resulfur);
            if(iich<0||iich>strlen(resulfur))continue;
            iimax=atchmax[iich];
            for(j=0;j<nres;j++)
            {
                if(i==j)continue;
                iiar=(int)(strchr(resaro,seq[j])-resaro);
                if(iiar>=0&&iiar<strlen(resaro))
                {
                    iimax1=atarmax[iiar];
                    aro=true;
                    if(seq[j]=='W')doub=true;
                    else doub=false;
                }
                else
                {
                    iiar=(int)(strchr(resbase,seq[j])-resbase);
                    if(iiar>=0&&iiar<strlen(resbase))
                    {
                        iimax1=atbamax[iiar];
                        aro=false;
                        if(seq[j]=='a'||seq[j]=='g'||seq[j]=='y')doub=true;
                        else doub=false;
                    }
                    else continue;
                }
                d=(float)NONEXIST;
                for(iip=0;iip<iimax;iip++)
                {
                    for(iiq=0;iiq<iimax1;iiq++)
                    {
                        dm=dist2(coordtot[i][iip],coordtot[j][iiq]);
                        if(dm<d)
                        {
                            d=dm;
                            iim=iip;
                            iima=iiq;
                        }
                    }
                }
                if(!doub)q=0;
                else
                {
                    d1=dist2(coordtot[i][iim],coord[j][0]);
                    d2=dist2(coordtot[i][iim],coord[j][1]);
                    if(d1<d2)q=0;
                    else q=1;
                }
                
                if(d<=(SPiDmax*SPiDmax))
                {
                    for(k=0;k<3;k++)ttN[k]=coordN[j][q][k];
                    if(seq[i]=='H')
                        for(k=0;k<3;k++)tt[k]=coord[i][0][k]-coord[j][q][k];
                    else
                        for(k=0;k<3;k++)tt[k]=coordtot[i][iim][k]-coord[j][q][k];
                    d1=dist2(tt,ttN);
                    for(k=0;k<3;k++)tt[k]=(-tt[k]);
                    d2=dist2(tt,ttN);
                    if(d2>d1)for(k=0;k<3;k++)tt[k]=(-tt[k]);
                    a=angle(tt,ttN,true);
                    
                    if(aro)anglemaxaro(d,&amax,coordtot[j][iima],coord[j][q],coordN[j][q],CATradcyl);
                    else anglemax(iiar,q,d,&amax,CATradcyl);
                    
                    if(a<=amax)
                    {
                        if(coordN[i][0][0]!=NONEXIST&&coordN[i][0][1]!=NONEXIST&&coordN[i][0][2]!=NONEXIST)
                        {
                            for(k=0;k<3;k++)tt[k]=coordN[i][0][k];
                            d1=dist2(tt,ttN);
                            for(k=0;k<3;k++)tt[k]=(-tt[k]);
                            d2=dist2(tt,ttN);
                            if(d2>d1)for(k=0;k<3;k++)tt[k]=(-tt[k]);
                            a1=angle(tt,ttN,true);
                            a1=a1*180.0/M_PI;
                        }
                        else a1=NONEXIST;
                        
                        strncpy(aux1,&hist[6*i],6);
                        strncpy(aux2,&hist[6*j],6);
                        
                        for(m=0;m<mem;m++)
                        {
                            if(seq[i]==memseq1[m]&&seq[j]==memseq2[m]&&
                               !strncmp(&aux1[1],&memhist1[m][1],5)&&
                               !strncmp(&aux2[1],&memhist2[m][1],5))break;
                        }
                        if(redundancy==0&&m<mem)continue;
                        if(flagADN[j]=='2')
                        {
                            for(m=0;m<mem;m++)
                            {
                                if(memflagADN2[m]=='2'&&seq[i]==memseq1[m]&&seq[j]==memseq2[m]&&
                                   !strncmp(&aux1[1],&memhist1[m][1],5)&&!strncmp(&memlongseq2[m*3],&longseq[j*3],3))break;
                            }
                            if(redundancy==0&&m<mem)continue;
                        }
                        
                        rep=false;
                        d=(float)sqrt((double)d);
                        a=a*180.0/M_PI;
                        amax=amax*180.0/M_PI;
                        fou=1;
                        
                        strcpy(memhist1[mem],aux1);
                        strcpy(memhist2[mem],aux2);
                        memseq1[mem]=seq[i];
                        memseq2[mem]=seq[j];
                        strncpy(&memlongseq2[mem*3],&longseq[j*3],3);
                        memprefix2[mem]=prefixes[j];
                        strncpy(auxlong,&longseq[j*3],3);
                        if(prefixes[j]=='\0')auxpref[0]=' ';
                        else auxpref[0]=prefixes[j];
                        memflagADN2[mem]=flagADN[j];
                        if(seq[i]!='K')mema1[mem]=a1;
                        repeat[mem]=rep;
                        
                        memdis[mem]=d;
                        memangle[mem]=a;
                        fl[mem]=fou;
                        
                        if(seq[j]=='H'||q==1)strcpy(mot,"C5");
                        else strcpy(mot,"C6");
                        
                        k=(int)(strchr(resulfur,seq[i])-resulfur);
                        strcpy(aux3,atsulfur[k][iim]);
                        if(flagADN[j]=='0')
                        {
                            k=(int)(strchr(resaro,seq[j])-resaro);
                            strcpy(aux4,ataro[k][iima]);
                        }
                        else strcpy(aux4,atbase[iima]);
                        
                        if(flagADN[j]=='2')
                        {
                            sprintf(ligne[mem],"  %c (%s) - %s %c [%s] (%s) : D=%4.2f\u212B {%3s-%3s} A=%4.1f° (%4.1f°) ",
                                    seq[i],aux1,auxlong,seq[j],mot,aux2,d,aux3,aux4,a,amax);
                            sprintf(lignecsv[mem],"%s,%c,%s,%s,%c,%s,%s,%.2f,%s,%s,%.1f,%.1f,",
                                    "sulfur-pi",seq[i],aux1,auxlong,seq[j],mot,aux2,d,aux3,aux4,a,amax);
                        }
                        else
                        {
                            sprintf(ligne[mem],"  %c (%s) - %c [%s] (%s) : D=%4.2f\u212B {%3s-%3s} A=%4.1f° (%4.1f°) ",
                                    seq[i],aux1,seq[j],mot,aux2,d,aux3,aux4,a,amax);
                            sprintf(lignecsv[mem],"%s,%c,%s,,%c,%s,%s,%.2f,%s,%s,%.1f,%.1f,",
                                    "sulfur-pi",seq[i],aux1,seq[j],mot,aux2,d,aux3,aux4,a,amax);
                        }
                        strcpy(memhb_seq[mem],hb_seq[i]);
                        strcpy(memhb_hist[mem],hb_hist[i]);
                        strcpy(memhb_at[mem],hb_at[i]);
                        for(k=0;k<strlen(hb_seq[i]);k++)
                        {
                            strncpy(aux3,&hb_hist[i][k*6],6);
                            strncpy(aux4,&hb_at[i][k*7],7);
                            if(k==0)sprintf(str," HB1: %c {%s} (%s)",hb_seq[i][k],aux4,aux3);
                            else sprintf(str,"; %c {%s} (%s)",hb_seq[i][k],aux4,aux3);
                            strcat(ligne[mem],str);
                            aux4[3]='\0';
                            sprintf(str,",%c,%s,%s,%s",hb_seq[i][k],aux3,&aux4[0],&aux4[4]);
                            strcat(lignecsv[mem],str);
                        }
                        if(strlen(hb_seq[i])>3)
                            printf("*** WARNING *** %c (%s) in %s makes more than 3 H-bonds\n",seq[i],aux1,prot);
                        if(flhbplus[0]==(char)true)
                        {
                            for(k=strlen(hb_seq[i]);k<3;k++)
                            {
                                sprintf(str,",,,,");
                                strcat(lignecsv[mem],str);
                            }
                        }
                        for(k=0;k<strlen(hb_seq[j]);k++)
                        {
                            strncpy(aux3,&hb_hist[j][k*6],6);
                            strncpy(aux4,&hb_at[j][k*7],7);
                            if(k==0)sprintf(str," HB2: %c {%s} (%s)",hb_seq[j][k],aux4,aux3);
                            else sprintf(str,"; %c {%s} (%s)",hb_seq[j][k],aux4,aux3);
                            strcat(ligne[mem],str);
                            aux4[3]='\0';
                            sprintf(str,",%c,%s,%s,%s",hb_seq[j][k],aux3,&aux4[0],&aux4[4]);
                            strcat(lignecsv[mem],str);
                        }
                        if(strlen(hb_seq[j])>3)
                            printf("*** WARNING *** %c (%s) in %s makes more than 3 H-bonds\n",seq[j],aux2,prot);
                        if(flhbplus[0]==(char)true)
                        {
                            for(k=strlen(hb_seq[j]);k<3;k++)
                            {
                                sprintf(str,",,,,");
                                strcat(lignecsv[mem],str);
                            }
                        }
                        strcat(ligne[mem],"\n");
                        strcat(lignecsv[mem],"\n");
                        mem++;
                        if(mem>=memcount)
                        {
                            memcountold=memcount;
                            memcount+=MEMMAX;
                            reallocate_mem(memcount,memcountold,&memseq1,&memseq2,&memlongseq2,&memflagADN2,
                                           &memprefix2,&memhist1,&memhist2,&ligne,&lignecsv,&memhb_hist,&memhb_seq,
                                           &memhb_at,&fl,&fldeja,&repeat,&mema1,&memdis,&memangle);
                        }
                    }
                }
            }
        }

        
        fprintf(seqout,"\n%s\n",prot);
                
        drap=(-1);
        for(m=memchargaro;m<memarosulfur;m++)
        {
            if(memflagADN2[m]!='0')continue;
            j=(int)(strchr(resaro,memseq2[m])-resaro);
            if(j<0||j>=strlen(resaro))continue;
            if(fl[m]<1)continue;
            if(drap<0)
            {
                fprintf(seqout," Protein-protein pi-pi interactions\n");
                drap=1;
            }
            fprintf(seqout,"%s",ligne[m]);
            fprintf(seqoutcsv,"%s,%s,%s",prot,"prot-prot",lignecsv[m]);
        }
    
        drap=(-1);
        for(m=0;m<memchargaro;m++)
        {
            if(memflagADN2[m]!='0')continue;
            j=(int)(strchr(resaro,memseq2[m])-resaro);
            if(j<0||j>=strlen(resaro))continue;
            
            if(fl[m]<1)continue;
            if(drap<0)
            {
                fprintf(seqout," Protein-protein cation-pi, His-pi and amino-pi interactions\n");
                drap=1;
            }
            fprintf(seqout,"%s",ligne[m]);
            fprintf(seqoutcsv,"%s,%s,%s",prot,"prot-prot",lignecsv[m]);
        }
    
        drap=(-1);
        for(m=memarosulfur;m<mem;m++)
        {
            if(memflagADN2[m]!='0')continue;
            j=(int)(strchr(resaro,memseq2[m])-resaro);
            if(j<0||j>=strlen(resaro))continue;
            
            if(fl[m]<1)continue;
            if(drap<0)
            {
                fprintf(seqout," Protein-protein sulfur-pi interactions\n");
                drap=1;
            }
            fprintf(seqout,"%s",ligne[m]);
            fprintf(seqoutcsv,"%s,%s,%s",prot,"prot-prot",lignecsv[m]);
        }
    
        drap=(-1);
        for(m=memchargaro;m<memarosulfur;m++)
        {
            if(memflagADN2[m]!='1')continue;
            j=(int)(strchr(resbase,memseq2[m])-resbase);
            if(j<0||j>=strlen(resbase))continue;
            
            if(fl[m]<1)continue;
            if(drap<0)
            {
                fprintf(seqout," Protein-DNA/RNA pi-pi interactions\n");
                drap=1;
            }
            fprintf(seqout,"%s",ligne[m]);
            fprintf(seqoutcsv,"%s,%s,%s",prot,"prot-DNA/RNA",lignecsv[m]);
        }
    
        drap=(-1);
        for(m=0;m<memchargaro;m++)
        {
            if(memflagADN2[m]!='1')continue;
            j=(int)(strchr(resbase,memseq2[m])-resbase);
            if(j<0||j>=strlen(resbase))continue;
            
            if(fl[m]<1)continue;
            if(drap<0)
            {
                fprintf(seqout," Protein-DNA/RNA cation-pi, His-pi and amino-pi interactions\n");
                drap=1;
            }
            fprintf(seqout,"%s",ligne[m]);
            fprintf(seqoutcsv,"%s,%s,%s",prot,"prot-DNA/RNA",lignecsv[m]);
        }
    
        drap=(-1);
        for(m=memarosulfur;m<mem;m++)
        {
            if(memflagADN2[m]!='1')continue;
            j=(int)(strchr(resbase,memseq2[m])-resbase);
            if(j<0||j>=strlen(resbase))continue;
            
            if(fl[m]<1)continue;
            if(drap<0)
            {
                fprintf(seqout," Protein-DNA/RNA sulfur-pi interactions\n");
                drap=1;
            }
            fprintf(seqout,"%s",ligne[m]);
            fprintf(seqoutcsv,"%s,%s,%s",prot,"prot-DNA/RNA",lignecsv[m]);
        }
    
        drap=(-1);
        for(m=memchargaro;m<memarosulfur;m++)
        {
            if(memflagADN2[m]!='2')continue;
            
            if(fl[m]<1)continue;
            if(drap<0)
            {
                fprintf(seqout," Protein-ligand pi-pi interactions\n");
                drap=1;
            }
            fprintf(seqout,"%s",ligne[m]);
            fprintf(seqoutcsv,"%s,%s,%s",prot,"prot-ligand",lignecsv[m]);
        }
    
        drap=(-1);
        for(m=0;m<memchargaro;m++)
        {
            if(memflagADN2[m]!='2')continue;

            if(fl[m]<1)continue;
            if(drap<0)
            {
                fprintf(seqout," Protein-ligand cation-pi, His-pi and amino-pi interactions\n");
                drap=1;
            }
            fprintf(seqout,"%s",ligne[m]);
            fprintf(seqoutcsv,"%s,%s,%s",prot,"prot-ligand",lignecsv[m]);
        }
    
        drap=(-1);
        for(m=memarosulfur;m<mem;m++)
        {
            if(memflagADN2[m]!='2')continue;

            if(fl[m]<1)continue;
            if(drap<0)
            {
                fprintf(seqout," Protein-ligand sulfur-pi interactions\n");
                drap=1;
            }
            fprintf(seqout,"%s",ligne[m]);
            fprintf(seqoutcsv,"%s,%s,%s",prot,"prot-ligand",lignecsv[m]);
        }
    
        comp=0;
        for(m=0;m<mem;m++)
        {
            if(memflagADN2[m]=='0')continue;
            for(k=0;k<strlen(memhb_seq[m]);k++)
            {
                sscanf(&memhb_hist[m][k*6+1],"%d",&i1);
                sscanf(&memhist2[m][1],"%d",&i2);
                if(memflagADN2[m]=='1'&&i1!=i2+1&&i2!=i1+1)continue;
                if(memflagADN2[m]=='2'&&i1!=i2)continue;
                if(memflagADN2[m]=='1'&&(memhb_at[m][k*7+6]=='*'||memhb_at[m][k*7+6]=='P'))
                    continue;
                i1=(int)(strchr(resbase,memhb_seq[m][k])-resbase);
                if(i1<0||i1>=strlen(resbase))continue;
                break;
            }
            if(k==strlen(memhb_seq[m]))continue;
            
            if(fl[m]<1)continue;
            if(comp==0)fprintf(seqout," Stair motifs\n");
            strncpy(auxh,&memhb_hist[m][k*6],6);
            strncpy(aux4,&memhb_at[m][k*7],7);
            fprintf(seqout,"  %c (%6s) - %c (%6s) - %c (%6s) {%s}\n",
                    memseq2[m],memhist2[m],memseq1[m],memhist1[m],memhb_seq[m][k],auxh,aux4);
            aux4[3]='\0';
            fprintf(seqoutcsv_chst,"%s,stair,%c,%s,%c,%s,%c,%s,%s,%s\n",
                    prot,memseq2[m],memhist2[m],memseq1[m],memhist1[m],memhb_seq[m][k],auxh,aux4,&aux4[4]);
            comp++;
        }
    
        deja=0;
        comp=0;
        for(m=0;m<mem;m++)
        {
            for(mm=m+1;mm<mem;mm++)
            {
                if(memseq1[m]==memseq1[mm]&&!strcmp(memhist1[m],memhist1[mm]))
                {
                    if(memseq2[m]==memseq2[mm]&&!strcmp(memhist2[mm],memhist2[m]))continue;
                    sprintf(ligne[deja],"  %c (%6s) - %c (%6s) - %c (%6s)",
                            memseq2[m],memhist2[m],memseq1[m],memhist1[m],memseq2[mm],memhist2[mm]);
                    sprintf(lignecsv[deja],"%c,%s,%c,%s,%c,%s",
                            memseq2[m],memhist2[m],memseq1[m],memhist1[m],memseq2[mm],memhist2[mm]);
                    if(fl[m]>=3&&fl[mm]>=3)fldeja[deja]=3;
                    else if(fl[m]>=2&&fl[mm]>=2)fldeja[deja]=2;
                    else fldeja[deja]=1;
                    for(t=0;t<deja;t++)
                        if(!strcmp(ligne[t],ligne[deja]))break;
                    
                    if((t==deja&&fldeja[deja]>=1)||
                       (t<deja&&fldeja[deja]>=1&&fldeja[t]<1))
                    {
                        if(comp==0)fprintf(seqout," Chains of pi interactions\n");
                        fprintf(seqout,"%s\n",ligne[deja]);
                        fprintf(seqoutcsv_chst,"%s,chain,%s\n",prot,lignecsv[deja]);
                        comp++;
                    }
                    deja++;
                }
                else if(memseq2[m]==memseq2[mm]&&!strcmp(memhist2[m],memhist2[mm]))
                {
                    if(memseq1[m]==memseq1[mm]&&!strcmp(memhist1[mm],memhist1[m]))continue;
                    sprintf(ligne[deja],"  %c (%6s) - %c (%6s) - %c (%6s)",
                            memseq1[m],memhist1[m],memseq2[m],memhist2[m],memseq1[mm],memhist1[mm]);
                    sprintf(lignecsv[deja],"%c,%s,%c,%s,%c,%s",
                            memseq1[m],memhist1[m],memseq2[m],memhist2[m],memseq1[mm],memhist1[mm]);
                    if(fl[m]>=3&&fl[mm]>=3)fldeja[deja]=3;
                    else if(fl[m]>=2&&fl[mm]>=2)fldeja[deja]=2;
                    else fldeja[deja]=1;
                    for(t=0;t<deja;t++)
                        if(!strcmp(ligne[t],ligne[deja]))break;
                    
                    if((t==deja&&fldeja[deja]>=1)||(t<deja&&fldeja[deja]>=1&&fldeja[t]<1))
                    {
                        if(comp==0)fprintf(seqout," Chains of pi interactions\n");
                        fprintf(seqout,"%s\n",ligne[deja]);
                        fprintf(seqoutcsv_chst,"%s,chain,%s\n",prot,lignecsv[deja]);
                        comp++;
                    }
                    deja++;
                }
                else if(memseq1[m]==memseq2[mm]&&!strcmp(memhist1[m],memhist2[mm]))
                {
                    if(memseq2[m]==memseq1[mm]&&!strcmp(memhist2[m],memhist1[mm]))
                        continue;
                    sprintf(ligne[deja],"  %c (%6s) - %c (%6s) - %c (%6s)",
                            memseq2[m],memhist2[m],memseq1[m],memhist1[m],memseq1[mm],memhist1[mm]);
                    sprintf(lignecsv[deja],"%c,%s,%c,%s,%c,%s",
                            memseq2[m],memhist2[m],memseq1[m],memhist1[m],memseq1[mm],memhist1[mm]);
                    if(fl[m]>=3&&fl[mm]>=3)fldeja[deja]=3;
                    else if(fl[m]>=2&&fl[mm]>=2)fldeja[deja]=2;
                    else fldeja[deja]=1;
                    for(t=0;t<deja;t++)
                        if(!strcmp(ligne[t],ligne[deja]))break;
                    
                    if((t==deja&&fldeja[deja]>=1)||
                       (t<deja&&fldeja[deja]>=1&&fldeja[t]<1))
                    {
                        if(comp==0)
                            fprintf(seqout," Chains of pi interactions\n");
                        fprintf(seqout,"%s\n",ligne[deja]);
                        fprintf(seqoutcsv_chst,"%s,chain,%s\n",prot,lignecsv[deja]);
                        comp++;
                    }
                    deja++;
                }
                else if(memseq2[m]==memseq1[mm]&&!strcmp(memhist2[m],memhist1[mm]))
                {
                    if(memseq1[m]==memseq2[mm]&&!strcmp(memhist1[m],memhist2[mm]))
                        continue;
                    sprintf(ligne[deja],"  %c (%6s) - %c (%6s) - %c (%6s)",
                            memseq1[m],memhist1[m],memseq2[m],memhist2[m],memseq2[mm],memhist2[mm]);
                    sprintf(lignecsv[deja],"%c,%s,%c,%s,%c,%s",
                            memseq1[m],memhist1[m],memseq2[m],memhist2[m],memseq2[mm],memhist2[mm]);
                    if(fl[m]>=3&&fl[mm]>=3)fldeja[deja]=3;
                    else if(fl[m]>=2&&fl[mm]>=2)fldeja[deja]=2;
                    else fldeja[deja]=1;
                    for(t=0;t<deja;t++)
                        if(!strcmp(ligne[t],ligne[deja]))break;
                    
                    if((t==deja&&fldeja[deja]>=1)||
                       (t<deja&&fldeja[deja]>=1&&fldeja[t]<1))
                    {
                        if(comp==0)
                            fprintf(seqout," Chains of pi interactions\n");
                        fprintf(seqout,"%s\n",ligne[deja]);
                        fprintf(seqoutcsv_chst,"%s,chain,%s\n",prot,lignecsv[deja]);
                        comp++;
                    }
                    deja++;
                }
            }
        }
        cc++;
    }

    printf("Found %d proteins\n",cc);
    sprintf(string,"rm %slist",repert);
    if(system(string)!=0)
    {
        printf("*** Command '%s' failed\n",string);
    }
        
    fprintf(seqout,"#---------------------------------------------------------------\n");
    fclose(seqout);
    fclose(seqoutcsv);
    fclose(seqlist);
}
/*---------------------------*/

