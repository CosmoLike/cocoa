
#include "clik.h"
#include "clik_helper.h"
#include <errno.h>
#include <stdio.h>
//#define _GNU_SOURCE
#include <string.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <poll.h>
#include <unistd.h>
#include <fcntl.h>
#include <limits.h>


void cfrd_pipe_err(int twinpip[2],error **err) {
  testErrorRetVA(pipe(twinpip) == -1,-98689,"Cannot create unnamed pipe (%s)",*err,__LINE__,,strerror(errno));
}
int cfrd_pollout_err(int fd, error **err) {
  struct pollfd fds;  
  int pollerr;
  
  fds.fd     = fd;
  fds.events = POLLOUT;
  pollerr  = poll(&fds,1,0);
  //_DEBUGHERE_("","");
  //_DEBUGHERE_("->%d",pollerr);
  testErrorRetVA(pollerr==0,-241342,"TIMEOUT ! >%.3f s",*err,__LINE__,pollerr,0);
  testErrorRetVA(pollerr<0,-241342,"Error on poll (%s)",*err,__LINE__,pollerr,strerror(errno));
  testErrorRetVA(fds.revents & POLLHUP,-241342,"Hangup ! (%d)",*err,__LINE__,pollerr,fds.revents);
  testErrorRetVA(fds.revents & POLLERR,-241342,"Poll err ! (%d)",*err,__LINE__,pollerr,fds.revents);
  return pollerr;
  
}
int cfrd_pollin_err(int fd, long delay, error **err) {
  struct pollfd fds;  
  int pollerr;
  
  if (delay<0) {
    return -1;
  }
  //_DEBUGHERE_("%d",delay);
  fds.fd     = fd;
  fds.events = POLLIN;
  pollerr  = poll(&fds,1,delay/1000);
  //_DEBUGHERE_("","");
  //_DEBUGHERE_("->%d",pollerr);
  testErrorRetVA(pollerr==0,-241342,"TIMEOUT ! >%.3f s",*err,__LINE__,pollerr,delay/1e6);
  testErrorRetVA(pollerr<0,-241342,"Error on poll (%s)",*err,__LINE__,pollerr,strerror(errno));
  testErrorRetVA(fds.revents & POLLHUP,-241342,"Hangup ! (%d)",*err,__LINE__,pollerr,fds.revents);
  testErrorRetVA(fds.revents & POLLERR,-241342,"Poll err ! (%d)",*err,__LINE__,pollerr,fds.revents);
  return pollerr;
  
}
void cfrd_send(int fd,char* mess,int nb,error **err) {
  cfrd_pollout_err(fd, err);
  forwardError(*err,__LINE__,);
  if (nb<0) {
    nb = strlen(mess);
  }
  //_DEBUGHERE_("SEND %d '%s' (%d %d)",fd,mess,nb,strlen(mess));
  testErrorRetVA(write(fd,mess,nb)==-1,-241342,"Error writing (%s)",*err,__LINE__,,strerror(errno));
}

int cfrd_read_err(int fd, char* buf,  error **err) {
  int res;
  int i;

  for(i=0;i<500;i++) {
    //cfrd_pollin_err(fd, delay, err);
    //forwardError(*err,__LINE__,);
    res = read(fd,&(buf[i]),1);
    testErrorRetVA(res==-1,-241342,"Error on read(%s)",*err,_LINE__,-1,strerror(errno));
    if (buf[i]=='\n') {
      buf[i] = '\0';
      return i;
    }     
  }
}

typedef struct {
  int pipin,pipout;
  } bicep;


double bicep_lkl(void* pbic, double* pars, error **err) {
  double lkl;
  bicep *bic;
  char cmd[500];
  bic = pbic;
  int i;

  for(i=0;i<601;i++) {
    sprintf(cmd,"%g\n",pars[i]);
    cfrd_send(bic->pipin,cmd,-1,err);
    forwardError(*err,__LINE__,0);
  }
  for(i=0;i<601;i++) {
    sprintf(cmd,"%g\n",pars[601+i]);
    cfrd_send(bic->pipin,cmd,-1,err);
    forwardError(*err,__LINE__,0);
  }
  while(1) {
    cfrd_read_err(bic->pipout,cmd,err);
    forwardError(*err,__LINE__,0);
    if (strcmp(cmd,"READY")==0) {
      break;
    }
  }
  cfrd_read_err(bic->pipout,cmd,err);
  forwardError(*err,__LINE__,0);
  
  sscanf(cmd,"%lg",&lkl);
  return lkl;

}

void free_bicep(void **pbic) {
  bicep *bic;
  error *_err,**err;

  _err = NULL;
  err = &_err;
  bic = *pbic;
  cfrd_send(bic->pipin,"stop\n",-1,err);
  forwardError(*err,__LINE__,);
  close(bic->pipin);
  close(bic->pipout);
}

cmblkl* clik_bicep_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  char directory_name[4096],pwd[4096],pwd2[4096],cmd[4096];
  cmblkl *cing;
  int pipin[2],pipout[2];
  int childpid;
  bicep * bic;

  bic = malloc_err(sizeof(bicep),err);
  forwardError(*err,__LINE__,NULL);
  
  cfrd_pipe_err(pipin,err); //pfd[0] is for reading, //pfd[1] is for writing
  forwardError(*err,__LINE__,NULL);
  cfrd_pipe_err(pipout,err); //pfd[0] is for reading, //pfd[1] is for writing
  forwardError(*err,__LINE__,NULL);

  childpid=fork();
  testErrorRet(childpid == -1,-241342,"fork failed",*err,__LINE__,NULL);
  
  /* this is the child */
  if (childpid == 0) {
    close(pipin[1]);
    close(pipout[0]);
    
    /* Duplicate stdin*/
    dup2(pipin[0],0);
    close(pipin[0]);   
    
    /* Duplicate stdout */
    dup2(pipout[1],1);
    close(pipout[1]);
    write(1,"test\n",5);
    execlp("bicep_main","bicep_main",(char*) NULL);
    // never returns from here, but we never know
    _DEBUGHERE_("bad oh bad","");
    exit(-1);
  } /* end of child process case */
  
  close(pipin[0]);
  close(pipout[1]);
  bic->pipin = pipin[1];
  bic->pipout = pipout[0];
  


  cldf_external(df,directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
 
  testErrorRetVA(getcwd(pwd2,4096)==NULL,-101010,"can't get cwd name (cause = '%s')",*err,__LINE__,NULL,strerror(errno));
  sprintf(cmd,"%s\n",pwd2);
  cfrd_send(bic->pipin,cmd,-1,err);
  forwardError(*err,__LINE__,NULL);
  
  cldf_external_cleanup(directory_name,pwd,err);  
  forwardError(*err,__LINE__,NULL);
  
  cing = init_cmblkl(bic, &bicep_lkl, 
                     &free_bicep,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);

  return cing;
}
