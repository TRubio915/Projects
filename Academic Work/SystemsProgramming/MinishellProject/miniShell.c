#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

void printCWD();
void parseArgs(char *buffer, char** args, int argsSize, int *nargs);
void exeProg(char** args, int argsSize);
int isIn(char** args, int argsSize, char* file);
int isOut(char** args, int argsSize, char* file);

int main()
{
  int bufSize = 20;
  char str[bufSize];
  char** argArr = (char**)malloc(20*sizeof(char*));
  int numFilled;
  while(1)
    {
      //Print Dir
      printCWD();
      //Get User input
      fgets(str, bufSize, stdin);
      //Parse input into tokens
      parseArgs(str,argArr,bufSize,&numFilled);
      //Exit
      if(strcmp(argArr[0],"exit") == 0)
	{
	  exit(0);
	}
      //Change Directory
      if(strcmp(argArr[0],"cd") == 0)
	{
	  if(chdir(argArr[1]) != 0)
	    {
	      printf("Directory Not Found\n");
	    }
	}
      //Executing a program
      else
	{
	  exeProg(argArr, numFilled);
	}
    }
}

//Function to print current working directory
void printCWD()
{
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  printf("\n%s> ", cwd);
}

// buffer is the input string of characters
// args is the output array of arguments.  It has already been created with argsSize slots.
// nargs as the number filled in that is passed back
void parseArgs(char *buffer, char** args, int argsSize, int *nargs)
{
  char *bufArgs[argsSize];
  char **cp;
  char *wbuf;
  int i, j;

  wbuf=buffer;
  bufArgs[0]=buffer;
  args[0]=buffer;

  for(cp=bufArgs; (*cp=strsep(&wbuf, " \n\t")) != NULL ;){
    if ((*cp != '\0') && (++cp >= &bufArgs[argsSize]))
      break;
  }

  for (j=i=0; bufArgs[i]!=NULL; i++){
    if(strlen(bufArgs[i])>0)
      args[j++]=bufArgs[i];
  }

  // Add the NULL as the end argument because we need that for later
  *nargs=j;
  args[j]=NULL;
}

void handler(int sig)
{
  printf("Previous Process Complete!\n");
}

//Function which utilizes forking to execute programs
void exeProg(char** args, int argsSize)
{
   
  int pid = fork();
  if(pid != 0)
    {
      if(strcmp(args[argsSize-1],"&") == 0)
	{
	  signal(SIGCHLD, handler);
	  args[argsSize-1] = "\0";
	  main();
	}
      else
	{
	  //parent
	  int status;
	  int result = waitpid(pid,&status,0); //Blocking to wait for child
	}
    }
  else
    {
      //child
      char* file = (char*)malloc(10*sizeof(char));//Filename to be input/output
      //If is output
      if(isOut(args,argsSize,file))
	{
	  freopen(file,"w",stdout);
	}
      //if is input
      else if(isIn(args,argsSize,file))
	{
	  freopen(file,"r",stdin);
	}
      //Run program
      int rc = execvp(args[0], args);
    }
}

//Check if input command
int isIn(char** args, int argsSize, char* file)
{
  for(int i = 0; i < argsSize; i++)
    {
      if(strcmp(args[i],"<") == 0)
	{
	  strcpy(file,args[i+1]);
	  for(int j = i; j < argsSize; j++)
	    {
	      args[j] = args[j+1];
	    }
	  argsSize--;
	  return 1;
	}
    }
  return 0;
}

//Check if output command
int isOut(char** args, int argsSize, char* file)
{
  for(int i = 0; i < argsSize; i++)
    {
      if(strcmp(args[i],">") == 0)
	{
	  strcpy(file,args[i+1]);
	  for(int j = i; j < argsSize; j++)
	    {
	      args[j] = args[j+1];
	    }
	  return 1;
	}
    }
  return 0;
}
