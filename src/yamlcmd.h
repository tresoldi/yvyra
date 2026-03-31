#ifndef YAMLCMD_H_
#define YAMLCMD_H_

/*
 *  yamlcmd.h: Convert parsed YAML config to yvyra NEXUS commands.
 */

/* Convert a YAML config file to a NEXUS string that can be executed.
   Returns a malloc'd string or NULL on error. Caller must free.
   Also populates charStateNames/charNStateNames globals if available. */
char *YamlToNexus (const char *yamlFilename);

/* Store per-character state names (called by YamlToNexus).
   charIdx is 0-based, stateNames is array of nStates strings. */
void SetCharStateNames (int charIdx, int nStates, const char **stateNames);

#endif /* YAMLCMD_H_ */
