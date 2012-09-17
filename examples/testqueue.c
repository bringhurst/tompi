#include <stdio.h>
#include <string.h>
#include "mpii.h"

int mymatch (MPII_Msg *want, MPII_Msg *have)
{
   if (want->type == have->type && want->req == have->req)
      return 1;
   else
      return 0;
}

int main (int argc, char *argv[])
{
   MPII_Msg_queue q;
   char s[512], *tok, cmd, type;
   int die, retry, ntok;
   MPII_Msg msg, lastsrch;

   MPII_queue_init (&q);

   while (1)
   {
      printf ("enqueue/search/repeat msg_avail/took_msg ptr\n");
      s[0] = '\0';
      fgets (s, 512, stdin);

      tok = strtok (s, " \t");
      if (tok == NULL)
         continue;
      cmd = tolower (*tok);
      ntok = 0;

      tok = strtok (NULL, " \t");
      if (tok == NULL)
         goto end_tok;
      type = tolower (*tok);
      die = 0;
      switch (type)
      {
         case 'm': msg.type = MSG_AVAIL; break;
         case 't': msg.type = TOOK_MSG; break;
         default:  printf ("Unknown type\n"); die = 1;
      }
      if (die)
         continue;
      ntok++;

      tok = strtok (NULL, " \t");
      if (tok == NULL)
         goto end_tok;
      msg.req = (void *) atoi (tok);
      ntok++;

      tok = strtok (NULL, " \t");
      if (tok != NULL)
      {
         printf ("Long line\n");
         continue;
      }

end_tok:
      switch (cmd)
      {
         case 'e':
            if (ntok != 2)
            {
               printf ("usage: enqueue type ptr\n");
               break;
            }
            MPII_enqueue (&q, &msg);
            break;
         case 's':
            retry = 0; /* fall through */
            if (ntok != 2)
            {
               printf ("usage: search type ptr\n");
               break;
            }
            lastsrch = msg;
         case 'r':
            MPII_queue_search (&retry, &q, mymatch, &lastsrch, &msg);
            break;
      }
   }
}

