#ifndef MATCH_PATH_SEPARATOR
#define MATCH_PATH_SEPARATOR


// Used with STL algorithms to find path separators
struct MatchPathSeparator
{
    bool operator()(char c) const
    {
        // Match both Unix-style and Windows-style path separators
        return c == '/' || c == '\\';
    }
};


#endif
