#ifndef YAMLPARSER_H_
#define YAMLPARSER_H_

/*
 *  yamlparser.h: Minimal YAML subset parser for yvyra configuration files.
 *
 *  Handles: scalars, block/flow mappings, block/flow sequences, comments,
 *  quoted strings, inline comments. No anchors, tags, or multi-document.
 */

#define YAML_MAX_DEPTH      32
#define YAML_MAX_CHILDREN   512
#define YAML_MAX_KEY_LEN    256
#define YAML_MAX_SCALAR_LEN 4096

typedef enum {
    YAML_SCALAR,
    YAML_MAPPING,
    YAML_SEQUENCE
} YamlNodeType;

typedef struct YamlNode {
    YamlNodeType    type;
    char            *key;           /* key for mapping entries (NULL for seq items) */
    char            *scalar;        /* value for SCALAR nodes */
    struct YamlNode **children;     /* child nodes for MAPPING/SEQUENCE */
    int             nChildren;
    int             maxChildren;
} YamlNode;

/* Parse a YAML file into a node tree. Returns root node or NULL on error. */
YamlNode *YamlParseFile (const char *filename);

/* Parse a YAML string into a node tree. Returns root node or NULL on error. */
YamlNode *YamlParseString (const char *text);

/* Free a node tree */
void YamlFreeNode (YamlNode *node);

/* Lookup helpers */
YamlNode *YamlGetChild (YamlNode *node, const char *key);
const char *YamlGetScalar (YamlNode *node, const char *key);
int YamlGetInt (YamlNode *node, const char *key, int defaultVal);
double YamlGetDouble (YamlNode *node, const char *key, double defaultVal);

/* Check if a scalar represents missing data (~) */
int YamlIsMissing (const char *scalar);

#endif /* YAMLPARSER_H_ */
