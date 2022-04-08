/************************************************************************/
/**

   Program:    sprotFTdist
   \file       sprotFTdist.c
   
   \version    V2.0
   \date       01.02.22
   \brief      SAAPdap plugin for swissprot feature distances
   
   \copyright  (c) UCL / Prof. Andrew C. R. Martin 2022
   \author     Prof. Andrew C. R. Martin
   \par
               Institute of Structural & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   \par
               andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified.

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0   01.02.22   Original   By: ACRM

*************************************************************************/
/* Debugging
*/
/* #define DEBUG         1 */
/* #define PRINTFEATURES 1 */


/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/MathType.h"


/************************************************************************/
/* Defines and macros
*/
#define MAXLABEL  8
#define SMALLBUFF 16
#define MAXBUFF   256
#define MAXSITE   100
#define BADCUTDIST (REAL)4.0

#define SERVERURL "http://www.bioinf.org.uk/servers/pdbsws/query.cgi?plain=1"
#define UNIPROTURL "https://www.uniprot.org/uniprot/"

#define PROGRESS(t) if(gVerbose) fprintf(stderr, "*** %s\n", t)

typedef struct _features
{
   REAL MinDistActSite,
        MinDistBinding,
        MinDistCABinding,
        MinDistDNABinding,
        MinDistNPBinding,
        MinDistMetal,
        MinDistModRes,
        MinDistCarbohyd,
        MinDistMotif,
        MinDistLipid;

   int  NActSite,
        NBinding,
        NCABinding,
        NDNABinding,
        NNPBinding,
        NMetal,
        NModRes,
        NCarbohyd,
        NMotif,
        NLipid;

   char ActSite[MAXSITE][MAXLABEL],
        Binding[MAXSITE][MAXLABEL],
        CABinding[MAXSITE][MAXLABEL],
        DNABinding[MAXSITE][MAXLABEL],
        NPBinding[MAXSITE][MAXLABEL],
        Metal[MAXSITE][MAXLABEL],
        ModRes[MAXSITE][MAXLABEL],
        Carbohyd[MAXSITE][MAXLABEL],
        Motif[MAXSITE][MAXLABEL],
        Lipid[MAXSITE][MAXLABEL];
}  FEATURES;

typedef struct _pdbsws
{
   char pdb[SMALLBUFF],
        chain[SMALLBUFF],
        resid[SMALLBUFF],
        pdbaa[SMALLBUFF],
        ac[SMALLBUFF],
        id[SMALLBUFF],
        upcount[SMALLBUFF],
        aa[SMALLBUFF];
   struct _pdbsws *next;
}  PDBSWS;


/************************************************************************/
/* Globals
*/
BOOL gVerbose = FALSE;


/************************************************************************/
/* Prototypes
*/
BOOL ParseCmdLine(int argc, char **argv, char *resid, char *newaa,
                  char *infile);
int main(int argc, char **argv);
void FindUniProtCode(char *pdbcode, char *resid, char *uniprotcode);
FEATURES FindFeatures(char *uniprotcode);
void MapFeaturesToPDB(FEATURES *features, char *upcode, char *pdbcode,
                      char *chain);
void MapFeature(char *label, char *upcode, char *pdbcode, char *chain,
                int nres, char resid[MAXSITE][MAXLABEL]);
void CalculateFeatureDistances(FEATURES *features, char *resid,
                               char *infile);
void Usage(void);
void CopyItem(char *body, char *key, char *dest);
PDBSWS *ParsePDBSWSResponse(char *response);
char *RunExternal(char *cmd);
void SetFeature(char *text, char *feature, int *nfeature,
                char residues[MAXSITE][MAXLABEL]);
int ExpandRange(char *range, int *residues);
void FindPDBResFromUniProt(char *upcode, char *upresid,
                           char *pdbcode, char *chain, char *resid);
void PopulateFeatureDistance(PDB *pdb, char *chain, int resnum,
                             char *insert, int nRes,
                             char resids[MAXSITE][MAXLABEL],
                             REAL *minDist);
void PrintAResult(char *label, REAL dist);
void PrintResults(FEATURES features);
void Info(void);
void ErrorExit(char *fmt, char *param);

#ifdef PRINTFEATURES
void PrintAFeature(char *label, int nres, char resids[MAXSITE][MAXLABEL],
                   REAL dist);
void PrintFeatures(FEATURES features);
#endif


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   \param[in]    argc  Argument count
   \param[in]    argv  Arguments
   \return             Status

   Main program

- 01.02.22  Original   By: ACRM
*/
int main(int argc, char **argv)
{
   char     resid[SMALLBUFF],
            newaa[SMALLBUFF],
            infile[MAXBUFF],
            uniprotcode[SMALLBUFF],
            pdbcode[SMALLBUFF];
   FEATURES features;
   
   
   if(ParseCmdLine(argc, argv, resid, newaa, infile))
   {
      char chain[MAXLABEL],
           insert[MAXLABEL];
      int  resnum;
      
      strcpy(pdbcode, blFNam2PDB(infile));
      PROGRESS("Finding UniProt code");
      FindUniProtCode(pdbcode, resid, uniprotcode);
#ifdef DEBUG
      printf("UP: %s\n", uniprotcode);
#endif
      PROGRESS("Finding features");
      features = FindFeatures(uniprotcode);
#ifdef PRINTFEATURES
      PrintFeatures(features);
#endif
      PROGRESS("Mapping features back to PDB");
      blParseResSpec(resid, chain, &resnum, insert);
      MapFeaturesToPDB(&features, uniprotcode, pdbcode, chain);
      PROGRESS("Calculating feature distances");
      CalculateFeatureDistances(&features, resid, infile);
#ifdef PRINTFEATURES
      PrintFeatures(features);
#endif
      PrintResults(features);
   }
   else
   {
      Usage();
   }
   
   return(0);
}


/************************************************************************/
/*>void MapFeaturesToPDB(FEATURES *features, char *upcode, char *pdbcode,
                         char *chain)
   ----------------------------------------------------------------------
*//**
   \param[in,out] *features Pointer to feature structure
   \param[in]     upcode    UniProt accession
   \param[in]     pdbcode   PDB code of interest 
   \param[in]     chain     PDB chain of interest

   Convert all features from UniProt number to PDB residue number (nnn[c])

- 01.02.22  Original   By: ACRM
*/
void MapFeaturesToPDB(FEATURES *features, char *upcode, char *pdbcode,
                      char *chain)
{
   char resid[SMALLBUFF];
   
   MapFeature("Active Site",  upcode, pdbcode, chain,
              features->NActSite,    features->ActSite);
   MapFeature("Binding",      upcode, pdbcode, chain,
              features->NBinding,    features->Binding);
   MapFeature("CA Binding",   upcode, pdbcode, chain,
              features->NCABinding,  features->CABinding);
   MapFeature("DNA Binding",  upcode, pdbcode, chain,
              features->NDNABinding, features->DNABinding);
   MapFeature("NP Binding",   upcode, pdbcode, chain,
              features->NNPBinding,  features->NPBinding);
   MapFeature("Metal",        upcode, pdbcode, chain,
              features->NMetal,      features->Metal);
   MapFeature("ModRes",       upcode, pdbcode, chain,
              features->NModRes,     features->ModRes);
   MapFeature("Carbohydrate", upcode, pdbcode, chain,
              features->NCarbohyd,   features->Carbohyd);
   MapFeature("Motif",        upcode, pdbcode, chain,
              features->NMotif,      features->Motif);
   MapFeature("Lipid",        upcode, pdbcode, chain,
              features->NLipid,      features->Lipid);
}


/************************************************************************/
/*>void MapFeature(char *label, char *upcode, char *pdbcode,
                   char *chain, int nres, char resid[MAXSITE][MAXLABEL])
   ---------------------------------------------------------
*//**
   \param[in]     label  A text label for debugging purposes
   \param[in]     upcode  UniProt accession
   \param[in]     pdbcode The PDB code of interest
   \param[in]     chain   The PDB chain of interest
   \param[in]     nres    The number of mapped residues in the feature
   \param[in,out] resid   Array of residue IDs for the mapped residues

   Convert an individual residue from UniProt number (nnn) to PDB number
   (nnn[c]) for the specified PDB code and chain

- 01.02.22  Original   By: ACRM
*/
void MapFeature(char *label, char *upcode, char *pdbcode,
                char *chain, int nres, char resid[MAXSITE][MAXLABEL])
{
   int i;

   for(i=0; i<nres; i++)
   {
      FindPDBResFromUniProt(upcode, resid[i], pdbcode, chain, resid[i]);
   }
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *resid, char *newaa,
                     char *infile)
   ------------------------------------------------------------------
*//**
   \param[in]    argc   Argument count
   \param[in]    argv   Arguments
   \param[out]   resid  The residue of interest
   \param[out]   newres The amino acid to which it is mutated
   \param[out]   infile The input PDB file
   \return              Success

   Parse the command line

- 01.02.22  Original   By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *resid, char *newaa,
                  char *infile)
{
   argc--;
   argv++;

   infile[0] = resid[0] = newaa[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         if(!strcmp(argv[0], "-vv"))
         {
            gVerbose = TRUE;
         }
         else if(!strcmp(argv[0], "-h"))
         {
            return(FALSE);
         }
         else if(!strcmp(argv[0], "-info"))
         {
            Info();
            return(TRUE);
         }
         else
         {
            return(FALSE);
         }
      }
      else
      {
         /* Check that there are only 3 arguments left                  */
         if(argc != 3)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(resid, argv[0]);
         /* Copy the first to infile                                    */
         strcpy(newaa, argv[1]);
         /* Copy the first to infile                                    */
         strcpy(infile, argv[2]);
         
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void FindUniProtCode(char *pdbcode, char *resid, char *uniprotcode)
   -------------------------------------------------------------------
*//**
   \param[in]    pdbcode      The PDB code
   \param[in]    resid        The residue identifier ([c]nnn[i])
   \param[out]   uniprotcode  The UniProt accession
   \return

   Find the UniProt code for a specified PDB code and residue ID

- 01.02.22  Original   By: ACRM
*/
void FindUniProtCode(char *pdbcode, char *resid, char *uniprotcode)
{
   char   chain[MAXLABEL],
          insert[MAXLABEL],
          url[MAXBUFF],
          *result,
          cmd[MAXBUFF];
   int    resnum;
   PDBSWS *pdbsws = NULL;
   
   blParseResSpec(resid, chain, &resnum, insert);
   sprintf(url, "%s&qtype=pdb&id=%s&chain=%s&res=%d%s",
           SERVERURL, pdbcode, chain, resnum, insert);
   KILLTRAILSPACES(url);
#ifdef DEBUG
   fprintf(stderr,"URL: %s\n", url);
#endif
   sprintf(cmd, "/usr/bin/curl -s '%s'", url);
   if((result = RunExternal(cmd))==NULL)
   {
      ErrorExit("External program failed: %s", cmd);
   }

   if((pdbsws = ParsePDBSWSResponse(result))!=NULL)
   {
      strcpy(uniprotcode, pdbsws->ac);
      FREELIST(pdbsws, PDBSWS);
      free(result);
   }
   else
   {
      free(result);
      ErrorExit("No data from PDBSWS call: %s", url);
   }
}


/************************************************************************/
/*>void FindPDBResFromUniProt(char *upcode, char *upresid, char *pdbcode,
                              char *chain, char *resid)
   ----------------------------------------------------------------------
*//**
   \param[in]   upcode   The UniProt code
   \param[in]   upresid  The residue number in the UniProt sequence
   \param[in]   pdbcode  The PDB code of interest
   \param[in]   chain    The PDB chain of interest
   \param[out]  resid    The residue number and insert code of the
                         PDB residue

   Find the PDB residue id (nnn[c]) for a specified UniProt code, UniProt
   residue number, PDB code and chain

- 01.02.22  Original   By: ACRM
*/
void FindPDBResFromUniProt(char *upcode, char *upresid, char *pdbcode,
                           char *chain, char *resid)
{
   char   url[MAXBUFF],
          cmd[MAXBUFF],
          *result;
   int    resnum;
   PDBSWS *pdbsws = NULL;
   
   sprintf(url, "%s&qtype=ac&id=%s&res=%s", SERVERURL, upcode, upresid);
   KILLTRAILSPACES(url);
#ifdef DEBUG
   fprintf(stderr,"URL: %s\n", url);
#endif
   sprintf(cmd, "/usr/bin/curl -s '%s'", url);
   if((result = RunExternal(cmd))==NULL)
   {
      if(gVerbose)
         fprintf(stderr, "No data from PDBSWS call: %s\n", url);
      return;
   }

   if((pdbsws=ParsePDBSWSResponse(result))!=NULL)
   {
      PDBSWS *p;
      for(p=pdbsws; p!=NULL; NEXT(p))
      {
         if(!strcmp(pdbcode, p->pdb) &&
            !strcmp(chain,   p->chain))
         {
            strcpy(resid, pdbsws->resid);
            break;
         }
      }
      FREELIST(pdbsws, PDBSWS);
   }
   else
   {
      if(gVerbose)
         fprintf(stderr, "No data from PDBSWS call: %s\n", url);
   }

   FREE(result);
}


/************************************************************************/
/*>FEATURES FindFeatures(char *uniprotcode)
   ----------------------------------------
*//**
   \param[in]   uniprotcode   A UniProt accession
   \return                    The relevant feature information

   Find the features from a UniProt entry

- 01.02.22  Original   By: ACRM
*/
FEATURES FindFeatures(char *uniprotcode)
{
   FEATURES features;
   struct http_response *httpResponse;
   char url[MAXBUFF];
   char cmd[MAXBUFF];
   char *result;

   sprintf(url, "%s%s.txt", UNIPROTURL, uniprotcode);
   sprintf(cmd, "/usr/bin/curl -s '%s'", url);

   if((result = RunExternal(cmd))==NULL)
   {
      ErrorExit("External program failed: %s", cmd);
   }

   features.NActSite    = 0;
   features.NBinding    = 0;
   features.NCABinding  = 0;
   features.NDNABinding = 0;
   features.NNPBinding  = 0;
   features.NMetal      = 0;
   features.NModRes     = 0;
   features.NCarbohyd   = 0;
   features.NMotif      = 0;
   features.NLipid      = 0;
   
   features.MinDistActSite    = (REAL)-1.0;
   features.MinDistBinding    = (REAL)-1.0;
   features.MinDistCABinding  = (REAL)-1.0;
   features.MinDistDNABinding = (REAL)-1.0;
   features.MinDistNPBinding  = (REAL)-1.0;
   features.MinDistMetal      = (REAL)-1.0;
   features.MinDistModRes     = (REAL)-1.0;
   features.MinDistCarbohyd   = (REAL)-1.0;
   features.MinDistMotif      = (REAL)-1.0;
   features.MinDistLipid      = (REAL)-1.0;
   
   SetFeature(result, "ACT_SITE", &(features.NActSite),
              features.ActSite);
   SetFeature(result, "BINDING",  &(features.NBinding),
              features.Binding);
   SetFeature(result, "CA_BIND",  &(features.NCABinding),
              features.CABinding);
   SetFeature(result, "DNA_BIND", &(features.NDNABinding),
              features.DNABinding);
   SetFeature(result, "NP_BIND",  &(features.NNPBinding),
              features.NPBinding);
   SetFeature(result, "METAL",    &(features.NMetal),
              features.Metal);
   SetFeature(result, "MOD_RES",  &(features.NModRes),
              features.ModRes);
   SetFeature(result, "CARBOHYD", &(features.NCarbohyd),
              features.Carbohyd);
   SetFeature(result, "MOTIF",    &(features.NMotif),
              features.Motif);
   SetFeature(result, "LIPID",    &(features.NLipid),
              features.Lipid);

   FREE(result);
   return(features);
}


/************************************************************************/
/*>void SetFeature(char *text, char *feature, int *nFtResidues,
                   char ftResidues[MAXSITE][MAXLABEL])
   ------------------------------------------------------------
*//**
   \param[in]   text         The whole UniProt record
   \param[in]   feature      The feature we are looking for
   \param[out]  *nFtResidues Number of residues for this feature
   \param[out]  ftResidues   Array of residue numbers in UniProt for
                             this feature

   Find an individual feature from a UniProt entry

- 01.02.22  Original   By: ACRM
*/
void SetFeature(char *text, char *feature, int *nFtResidues,
                char ftResidues[MAXSITE][MAXLABEL])
{
   static char *buffer = NULL;
   char *cstart,
        *cstop,
        key[MAXBUFF];

   if(buffer==NULL)
   {
      if((buffer = (char *)malloc((strlen(text)+2)*sizeof(char)))==NULL)
      {
         ErrorExit("No memory for buffer", NULL);
      }
   }
   strcpy(buffer, text);

   sprintf(key, "FT   %s", feature);

   cstart=cstop=buffer;
   while((cstart!=NULL) && (cstop!=NULL) && (*cstart != '\0'))
   {
      if((cstop = strchr(cstart, '\n'))!=NULL)
      {
         *cstop = '\0';
      }

      if(!strncmp(cstart, key, strlen(key)))
      {
         int residues[MAXSITE];
         int nInRange;
         int i;
            
         nInRange = ExpandRange(cstart, residues);
         for(i=0; i<nInRange; i++)
         {
            char resid[MAXLABEL];
            sprintf(resid, "%d", residues[i]);
            strcpy(ftResidues[*nFtResidues], resid);
            (*nFtResidues)++;
         }
      }
      if(cstop!=NULL)
         *cstop = '\n';
      
      cstart = cstop+1;
   }
}


/************************************************************************/
/*>int ExpandRange(char *range, int *residues)
   -------------------------------------------
*//**
   \param[in]   range    Residue number or range (xxx..xxx)
   \param[out]  residues Array of residue numbers
   \return               Number of residues in range

   Takes a string containing an number, or a range of numbers and 
   fills in an int array with all the numbers in the range

- 01.02.22  Original   By: ACRM
*/
int ExpandRange(char *range, int *residues)
{
   int  nResidues = 0,
        start     = 0,
        stop      = 0,
        i;
   char *dotdot,
        *value;
   
   value = strrchr(range, ' ') + 1;
   
   if((dotdot=strstr(value, ".."))!=NULL)
   {
      *dotdot = '\0';
      start   = atoi(value);
      stop    = atoi(dotdot+2);
      for(i=start; i<=stop; i++)
      {
         residues[nResidues++] = i;
      }
   }
   else
   {
      residues[0] = atoi(value);
      nResidues   = 1;
   }

   return(nResidues);
}


/************************************************************************/
/*>void CalculateFeatureDistances(FEATURES *features, char *resid,
                                  char *infile)
   ---------------------------------------------------------------
*//**
   \param[in,out]  features    Feature information
   \param[in]      resid       Key residue of interest
   \param[in]      infile      PDB file

   Fill in the FEATURES structure with the closest distance from a
   key residue to all residues in each feature type

- 01.02.22  Original   By: ACRM
*/
void CalculateFeatureDistances(FEATURES *features, char *resid,
                               char *infile)
{
   FILE *fp  = NULL;
   PDB  *pdb = NULL;
   char chain[MAXLABEL],
        insert[MAXLABEL];
   int  resnum,
        natoms;

   blParseResSpec(resid, chain, &resnum, insert);
   
   if((fp=fopen(infile, "r"))!=NULL)
   {
      if((pdb=blReadPDBAtoms(fp, &natoms))!=NULL)
      {
         PopulateFeatureDistance(pdb, chain, resnum, insert,
                                 features->NActSite,
                                 features->ActSite,
                                 &(features->MinDistActSite));
         PopulateFeatureDistance(pdb, chain, resnum, insert,
                                 features->NBinding,
                                 features->Binding,
                                 &(features->MinDistBinding));
         PopulateFeatureDistance(pdb, chain, resnum, insert,
                                 features->NCABinding,
                                 features->CABinding,
                                 &(features->MinDistCABinding));
         PopulateFeatureDistance(pdb, chain, resnum, insert,
                                 features->NDNABinding,
                                 features->DNABinding,
                                 &(features->MinDistDNABinding));
         PopulateFeatureDistance(pdb, chain, resnum, insert,
                                 features->NNPBinding,
                                 features->NPBinding,
                                 &(features->MinDistNPBinding));
         PopulateFeatureDistance(pdb, chain, resnum, insert,
                                 features->NMetal,
                                 features->Metal,
                                 &(features->MinDistMetal));
         PopulateFeatureDistance(pdb, chain, resnum, insert,
                                 features->NModRes,
                                 features->ModRes,
                                 &(features->MinDistModRes));
         PopulateFeatureDistance(pdb, chain, resnum, insert,
                                 features->NCarbohyd,
                                 features->Carbohyd,
                                 &(features->MinDistCarbohyd));
         PopulateFeatureDistance(pdb, chain, resnum, insert,
                                 features->NMotif,
                                 features->Motif,
                                 &(features->MinDistMotif));
         PopulateFeatureDistance(pdb, chain, resnum, insert,
                                 features->NLipid,
                                 features->Lipid,
                                 &(features->MinDistLipid));
      }
      else
      {
         ErrorExit("No atoms read from PDB file: %s", infile);
      }
   }
   else
   {
      ErrorExit("Unable to open file: %s", infile);
   }
}


/************************************************************************/
/*>PDBSWS *ParsePDBSWSResponse(char *response)
   -------------------------------------------
*//**
   \param[in]   response   Results from a PDBSWS query
   \return                 Linked list of PDBSWS records

   Take the string response from a PDBSWS query and convert it into
   a parsed linked list of PDBSWS structures

- 01.02.22  Original   By: ACRM
*/
PDBSWS *ParsePDBSWSResponse(char *response)
{
   char   *slashslash = NULL;
   PDBSWS *pdbsws     = NULL,
          *p          = NULL;
   BOOL   more        = FALSE;

   if(response == NULL)
      return(NULL);
   
   if(strstr(response, "400 Bad Request"))
   {
      ErrorExit("400 Bad Request", NULL);
   }
   
   do{
      if(p==NULL)
      {
         INIT(pdbsws, PDBSWS);
         p = pdbsws;
      }
      else
      {
         ALLOCNEXT(p, PDBSWS);
      }
      if(p==NULL)
      {
         ErrorExit("No memory for PDBSWS data", NULL);
      }
      
      CopyItem(response, "PDB: ",     p->pdb);
      CopyItem(response, "CHAIN: ",   p->chain);
      CopyItem(response, "RESID: ",   p->resid);
      CopyItem(response, "PDBAA: ",   p->pdbaa);
      CopyItem(response, "AC: ",      p->ac);
      CopyItem(response, "ID: ",      p->id);
      CopyItem(response, "UPCOUNT: ", p->upcount);
      CopyItem(response, "AA: ",      p->aa);

      /* Move on to next entry                                          */
      more = FALSE;
      if((slashslash = strstr(response, "//"))!=NULL)
         response = slashslash+2;

      if((slashslash = strstr(response, "//"))!=NULL)
         more = TRUE;
   }  while(more);

   return(pdbsws);
}


/************************************************************************/
/*>void CopyItem(char *response, char *key, char *dest)
   ----------------------------------------------------
*//**
   \param[in]    response   Complete UniProt record
   \param[in]    key        Record type of interest
   \param[out]   dest       Destination for data associated with key

   Extract the data for a given key (e.g. 'FT   BINDING') and extracts
   the associated residue number information

- 01.02.22  Original   By: ACRM
*/
void CopyItem(char *response, char *key, char *dest)
{
   char buffer[MAXBUFF],
        *chp;

   if(dest!=NULL)
   {
      *dest = '\0';
      if((chp=strstr(response, key))!=NULL)
      {
         chp += strlen(key);
         
         strncpy(buffer, chp, MAXBUFF-1);
         if((chp = strchr(buffer, '\n'))!=NULL)
         {
            *chp = '\0';
         }
         
         strncpy(dest, buffer, SMALLBUFF-1);
      }
   }
}


/************************************************************************/
/*>char *RunExternal(char *cmd)
   ----------------------------
*//**
   \param[in]    cmd   Command to be run
   \return             Results (malloc'd)

   Runs an external program returning the results in a malloc'd string

- 01.02.22  Original   By: ACRM
*/
char *RunExternal(char *cmd)
{
   FILE *fp;
   char buffer[MAXBUFF],
        *result = NULL;
   
   if ((fp=popen(cmd, "r")) == NULL)
   {
      return(NULL);
   }
   
   while (fgets(buffer, MAXBUFF, fp) != NULL)
   {
      result=blStrcatalloc(result, buffer);
   }

   if(pclose(fp))
   {
      return(NULL);
   }
   
   return(result);
}


/************************************************************************/
/*>void PopulateFeatureDistance(PDB *pdb,
                                char *chain, int resnum, char *insert,
                                int nRes, char resids[MAXSITE][MAXLABEL],
                                REAL *minDist)
   ----------------------------------------------------------------------
*//**
   \param[in]    pdb     PDB linked list
   \param[in]    chain   PDB chain of interest
   \param[in]    resnum  PDB resnum of key residue
   \param[in]    insert  PDB insert of key residue
   \param[in]    nRes    Number of feature residues
   \param[in]    resids  Res IDs (nnn[c]) of feature residues
   \param[out]   minDist Minimum distance to a feature residue (-1 if no
                         feature residues for this feature)


   Calculates the minimum distance between a key residue and a set of
   other residues.

- 01.02.22  Original   By: ACRM
*/
void PopulateFeatureDistance(PDB *pdb,
                             char *chain, int resnum, char *insert,
                             int nRes, char resids[MAXSITE][MAXLABEL],
                             REAL *minDist)
{
   PDB *keyRes, *keyResNext,
       *ftRes, *ftResNext;

   if(nRes)
   {
      if((keyRes = blFindResidue(pdb, chain, resnum, insert))!=NULL)
      {
         REAL minFtDist = (REAL)100000.0;
         int i;
         keyResNext = blFindNextResidue(keyRes);

         for(i=0; i<nRes; i++)
         {
            char ftChain[MAXLABEL], ftInsert[MAXLABEL];
            int  ftResnum;
            blParseResSpec(resids[i], ftChain, &ftResnum, ftInsert);
            if((ftRes = blFindResidue(pdb, chain, ftResnum, ftInsert))
               != NULL)
            {
               PDB  *p, *q;
               
               ftResNext = blFindNextResidue(ftRes);
               for(p=keyRes; p!=keyResNext; NEXT(p))
               {
                  for(q=ftRes; q!=ftResNext; NEXT(q))
                  {
                     REAL dist = DIST(p, q);
                     if(dist < minFtDist)
                     {
                        minFtDist = dist;
                     }
                  }
               }
            }
         }

         if(((*minDist < 0.0) || (minFtDist < *minDist)) &&
            (minFtDist < (99999.0)))
         {
            *minDist = minFtDist;
         }
      }
   }
}


/************************************************************************/
/*>void PrintResults(FEATURES features)
   ------------------------------------
*//**
   \param[in]   features   Feature information

   Print the results as a JSON string

- 01.02.22  Original   By: ACRM
*/
void PrintResults(FEATURES features)
{
   if(((features.MinDistActSite > (-0.5))    &&
       (features.MinDistActSite < BADCUTDIST))    ||
      ((features.MinDistBinding > (-0.5))    &&
       (features.MinDistBinding < BADCUTDIST))    ||
      ((features.MinDistCABinding > (-0.5))  &&
       (features.MinDistCABinding < BADCUTDIST))  ||
      ((features.MinDistDNABinding > (-0.5)) &&
       (features.MinDistDNABinding < BADCUTDIST)) ||
      ((features.MinDistNPBinding > (-0.5))  &&
       (features.MinDistNPBinding < BADCUTDIST))  ||
      ((features.MinDistMetal > (-0.5))      &&
       (features.MinDistMetal < BADCUTDIST))      ||
      ((features.MinDistModRes > (-0.5))     &&
       (features.MinDistModRes < BADCUTDIST))     ||
      ((features.MinDistCarbohyd > (-0.5))   &&
       (features.MinDistCarbohyd < BADCUTDIST))   ||
      ((features.MinDistMotif > (-0.5))      &&
       (features.MinDistMotif < BADCUTDIST))      ||
      ((features.MinDistLipid > (-0.5))      &&
       (features.MinDistLipid < BADCUTDIST)))
   {
      printf("{\"SprotFTdist-BOOL\": \"BAD\"");
   }
   else
   {
      printf("{\"SprotFTdist-BOOL\": \"OK\"");
   }
      
   PrintAResult("SprotFTdist-ACT_SITE", features.MinDistActSite);
   PrintAResult("SprotFTdist-BINDING",  features.MinDistBinding);
   PrintAResult("SprotFTdist-CA_BIND",  features.MinDistCABinding);
   PrintAResult("SprotFTdist-DNA_BIND", features.MinDistDNABinding);
   PrintAResult("SprotFTdist-NP_BIND",  features.MinDistNPBinding);
   PrintAResult("SprotFTdist-METAL",    features.MinDistMetal);
   PrintAResult("SprotFTdist-MOD_RES",  features.MinDistModRes);
   PrintAResult("SprotFTdist-CARBOHYD", features.MinDistCarbohyd);
   PrintAResult("SprotFTdist-MOTIF",    features.MinDistMotif);
   PrintAResult("SprotFTdist-LIPID",    features.MinDistLipid);

   printf("}\n");
}


/************************************************************************/
/*>void PrintAResult(char *label, REAL dist)
   -----------------------------------------
*//**
   \param[in]   label   Descriptor of this feature result
   \param[in]   dist    Minimum distance to this feature

   Print an individual result in JSON format

- 01.02.22  Original   By: ACRM
*/
void PrintAResult(char *label, REAL dist)
{
   int i;
   
   printf(", \"%s\": \"%.3f\"", label, dist);
}


#ifdef PRINTFEATURES
/************************************************************************/
/*>void PrintFeatures(FEATURES features)
   -------------------------------------
*//**
   \param[in]   features   Feature information

   Print the feature information for debugging purposes

- 01.02.22  Original   By: ACRM
*/
void PrintFeatures(FEATURES features)
{
   PrintAFeature("SprotFTdist-ACT_SITE", features.NActSite,
                 features.ActSite,       features.MinDistActSite);
   PrintAFeature("SprotFTdist-BINDING",  features.NBinding,
                 features.Binding,       features.MinDistBinding);
   PrintAFeature("SprotFTdist-CA_BIND",  features.NCABinding,
                 features.CABinding,     features.MinDistCABinding);
   PrintAFeature("SprotFTdist-DNA_BIND", features.NDNABinding,
                 features.DNABinding,    features.MinDistDNABinding);
   PrintAFeature("SprotFTdist-NP_BIND",  features.NNPBinding,
                 features.NPBinding,     features.MinDistNPBinding);
   PrintAFeature("SprotFTdist-METAL",    features.NMetal,
                 features.Metal,         features.MinDistMetal);
   PrintAFeature("SprotFTdist-MOD_RES",  features.NModRes,
                 features.ModRes,        features.MinDistModRes);
   PrintAFeature("SprotFTdist-CARBOHYD", features.NCarbohyd,
                 features.Carbohyd,      features.MinDistCarbohyd);
   PrintAFeature("SprotFTdist-MOTIF",    features.NMotif,
                 features.Motif,         features.MinDistMotif);
   PrintAFeature("SprotFTdist-LIPID",    features.NLipid,
                 features.Lipid,         features.MinDistLipid);
}


/************************************************************************/
/*>void PrintAFeature(char *label, int nres, 
                      char resids[MAXSITE][MAXLABEL], REAL dist)
   -------------------------------------------------------------
*//**
   \param[in]   label    Label for this feature
   \param[in]   nres     Number of residues for this feature
   \param[in]   resids   Residue IDs for this feature
   \param[in]   dist     Minimum distance to this feature

   Print an individual feature for debugging purposes

- 01.02.22  Original   By: ACRM
*/
void PrintAFeature(char *label, int nres, char resids[MAXSITE][MAXLABEL],
                   REAL dist)
{
   int i;
   
   fprintf(stderr, "%s: (%d)", label, nres);
   for(i=0; i<nres; i++)
   {
      fprintf(stderr, " %s", resids[i]);
   }
   fprintf(stderr, " [%.3f]\n", dist);
}
#endif


/************************************************************************/
/*>void Info(void)
   ---------------
*//**

   Prints the information string for the plugin

- 01.02.22  Original   By: ACRM
*/
void Info(void)
{
   printf("Finding distances to SwissProt features\n");
   exit(0);
}


/************************************************************************/
/*>void ErrorExit(char *fmt, char *param)
   --------------------------------------
*//**
   \param[in]    fmt    String to print (may have a %s)
   \param[in]    param  String to insert in fmt (or NULL)

   Prints an error message in JSON format and exits

- 01.02.22  Original   By: ACRM
*/
void ErrorExit(char *fmt, char *param)
{
   printf("{\"SprotFTdist-ERROR\": \"");
   if(param != NULL)
   {
      printf(fmt, param);
   }
   else
   {
      printf(fmt);
   }
   printf("\"}\n");
   exit(1);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

   Prints a usage message

- 01.02.22  Original   By: ACRM
*/
void Usage(void)
{
   printf("\nsprotFTdist V2.0 (c) UCL, Prof. Andrew C.R. Martin, \
Barbara A. Mikucka\n");

   printf("\nUsage: sprotfeatures.py [-vv][-nocache][-force][-info]\n");
   printf("       [chain]resnum[insert] newaa pdbfile\n");

   printf("\n       (newaa maybe 3-letter or 1-letter code)\n");

   printf("\n       -vv      Verbose\n");
   printf("       -nocache Do not cache results\n");
   printf("       -force   Force calculation even if results are \
cached\n");
   printf("       -info    Prints a 1-line summary of what the plugin \
is doing\n");
   printf("                and exits\n");

   printf("\nCalculates the distances of a mutant residue to the \
closest of each\n");
   printf("SwissProt feature type.\n\n");
}
