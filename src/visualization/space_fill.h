#ifndef _SPACE_FILL_H_
#define _SPACE_FILL_H_

#include "../utility/Types.h"
#include "../utility/Vector3d.h"

class space_fill{
private:
    std::vector<CVector3d> m_atoms;
    std::vector<double> m_radii;
    std::vector<AtomType> m_types;
public:
    space_fill();
    space_fill(std::vector<CVector3d> atoms, std::vector<double> radii, std::vector<AtomType> types);
    ~space_fill();

    void write();

    void write_all_atoms();
    void write_atoms(AtomType at, std::string sfx);
};

// sphere data
const std::vector<double> points = {  
    0, 0, -1, 
    0.203181, -0.147618, -0.96795, 
    -0.077607, -0.238853, -0.96795, 
    0.723607, -0.525725, -0.44722, 
    0.609547, -0.442856, -0.657519, 
    0.812729, -0.295238, -0.502301, 
    -0.251147, 0, -0.967949, 
    -0.077607, 0.238853, -0.96795, 
    0.203181, 0.147618, -0.96795, 
    0.860698, -0.442858, -0.251151, 
    -0.276388, -0.850649, -0.44722, 
    -0.029639, -0.864184, -0.502302, 
    -0.155215, -0.955422, -0.251152, 
    -0.894426, 0, -0.447216, 
    -0.831051, -0.238853, -0.502299, 
    -0.956626, -0.147618, -0.251149, 
    -0.276388, 0.850649, -0.44722, 
    -0.483971, 0.716565, -0.502302, 
    -0.436007, 0.864188, -0.251152, 
    0.723607, 0.525725, -0.44722, 
    0.531941, 0.681712, -0.502302, 
    0.687159, 0.681715, -0.251152, 
    0.687159, -0.681715, -0.251152, 
    -0.436007, -0.864188, -0.251152, 
    -0.956626, 0.147618, -0.251149, 
    -0.155215, 0.955422, -0.251152, 
    0.860698, 0.442858, -0.251151, 
    0.276388, -0.850649, 0.44722, 
    0.483971, -0.716565, 0.502302, 
    0.232822, -0.716563, 0.657519, 
    -0.723607, -0.525725, 0.44722, 
    -0.531941, -0.681712, 0.502302, 
    -0.609547, -0.442856, 0.657519, 
    -0.723607, 0.525725, 0.44722, 
    -0.812729, 0.295238, 0.502301, 
    -0.609547, 0.442856, 0.657519, 
    0.276388, 0.850649, 0.44722, 
    0.029639, 0.864184, 0.502302, 
    0.232822, 0.716563, 0.657519, 
    0.894426, 0, 0.447216, 
    0.831051, 0.238853, 0.502299, 
    0.753442, 0, 0.657515, 
    0.251147, 0, 0.967949, 
    0.077607, 0.238853, 0.96795, 
    0, 0, 1, 
    0.52573, 0, 0.850652, 
    0.3618, 0.262863, 0.894429, 
    0.638194, 0.262864, 0.72361, 
    0.162456, 0.499995, 0.850654, 
    0.447209, 0.525728, 0.723612, 
    0.688189, 0.499997, 0.525736, 
    0.483971, 0.716565, 0.502302, 
    -0.203181, 0.147618, 0.96795, 
    -0.138197, 0.425319, 0.89443, 
    -0.05279, 0.688185, 0.723612, 
    -0.425323, 0.309011, 0.850654, 
    -0.361804, 0.587778, 0.723612, 
    -0.262869, 0.809012, 0.525738, 
    -0.531941, 0.681712, 0.502302, 
    -0.203181, -0.147618, 0.96795, 
    -0.44721, 0, 0.894429, 
    -0.670817, 0.162457, 0.723611, 
    -0.425323, -0.309011, 0.850654, 
    -0.670817, -0.162457, 0.723611, 
    -0.850648, 0, 0.525736, 
    -0.812729, -0.295238, 0.502301, 
    0.077607, -0.238853, 0.96795, 
    -0.138197, -0.425319, 0.89443, 
    -0.361804, -0.587778, 0.723612, 
    0.162456, -0.499995, 0.850654, 
    -0.05279, -0.688185, 0.723612, 
    -0.262869, -0.809012, 0.525738, 
    0.029639, -0.864184, 0.502302, 
    0.3618, -0.262863, 0.894429, 
    0.447209, -0.525728, 0.723612, 
    0.638194, -0.262864, 0.72361, 
    0.688189, -0.499997, 0.525736, 
    0.831051, -0.238853, 0.502299, 
    0.956626, 0.147618, 0.251149, 
    0.951058, 0.309013, 0, 
    0.861804, 0.425322, 0.276396, 
    0.809019, 0.587782, 0, 
    0.670821, 0.688189, 0.276397, 
    0.587786, 0.809017, -0, 
    0.436007, 0.864188, 0.251152, 
    0.155215, 0.955422, 0.251152, 
    0, 1, 0, 
    -0.138199, 0.951055, 0.276397, 
    -0.309016, 0.951057, -0, 
    -0.447215, 0.850649, 0.276397, 
    -0.587786, 0.809017, -0, 
    -0.687159, 0.681715, 0.251152, 
    -0.860698, 0.442858, 0.251151, 
    -0.951058, 0.309013, 0, 
    -0.947213, 0.162458, 0.276396, 
    -1, -0, 1e-06, 
    -0.947213, -0.162458, 0.276397, 
    -0.951058, -0.309013, -0, 
    -0.860698, -0.442858, 0.251151, 
    -0.687159, -0.681715, 0.251152, 
    -0.587786, -0.809017, 0, 
    -0.447216, -0.850648, 0.276397, 
    -0.309017, -0.951056, -1e-06, 
    -0.138199, -0.951055, 0.276397, 
    0, -1, -0, 
    0.155215, -0.955422, 0.251152, 
    0.436007, -0.864188, 0.251152, 
    0.587786, -0.809017, 0, 
    0.67082, -0.68819, 0.276396, 
    0.809019, -0.587783, -2e-06, 
    0.861804, -0.425323, 0.276394, 
    0.951058, -0.309013, -0, 
    0.956626, -0.147618, 0.251149, 
    0.309017, 0.951056, -0, 
    0.447216, 0.850648, -0.276398, 
    0.138199, 0.951055, -0.276398, 
    0.262869, 0.809012, -0.525738, 
    -0.029639, 0.864184, -0.502302, 
    -0.809018, 0.587783, -0, 
    -0.670819, 0.688191, -0.276397, 
    -0.861803, 0.425324, -0.276396, 
    -0.688189, 0.499997, -0.525736, 
    -0.831051, 0.238853, -0.502299, 
    -0.809018, -0.587783, 0, 
    -0.861803, -0.425324, -0.276396, 
    -0.670819, -0.688191, -0.276397, 
    -0.688189, -0.499997, -0.525736, 
    -0.483971, -0.716565, -0.502302, 
    0.309017, -0.951056, 0, 
    0.138199, -0.951055, -0.276398, 
    0.447216, -0.850648, -0.276398, 
    0.262869, -0.809012, -0.525738, 
    0.531941, -0.681712, -0.502302, 
    1, 0, 0, 
    0.947213, -0.162458, -0.276396, 
    0.947213, 0.162458, -0.276396, 
    0.850648, 0, -0.525736, 
    0.812729, 0.295238, -0.502301, 
    0.609547, 0.442856, -0.657519, 
    0.425323, 0.309011, -0.850654, 
    0.361803, 0.587779, -0.723612, 
    0.138197, 0.425321, -0.894429, 
    0.052789, 0.688186, -0.723611, 
    -0.162456, 0.499995, -0.850654, 
    -0.232822, 0.716563, -0.657519, 
    -0.447211, 0.525727, -0.723612, 
    -0.361801, 0.262863, -0.894429, 
    -0.638195, 0.262863, -0.723609, 
    -0.52573, 0, -0.850652, 
    -0.753442, 0, -0.657515, 
    -0.638195, -0.262864, -0.723609, 
    -0.361801, -0.262864, -0.894428, 
    -0.447211, -0.525729, -0.72361, 
    -0.162456, -0.499995, -0.850654, 
    -0.232822, -0.716563, -0.657519, 
    0.670817, 0.162457, -0.723611, 
    0.670818, -0.162458, -0.72361, 
    0.447211, -1e-06, -0.894428, 
    0.425323, -0.309011, -0.850654, 
    0.05279, -0.688185, -0.723612, 
    0.138199, -0.425321, -0.894429, 
    0.361805, -0.587779, -0.723611
};

const std::vector<int> faces = {
    0, 1, 2, 
    3, 4, 5, 
    0, 2, 6, 
    0, 6, 7, 
    0, 7, 8, 
    3, 5, 9, 
    10, 11, 12, 
    13, 14, 15, 
    16, 17, 18, 
    19, 20, 21, 
    3, 9, 22, 
    10, 12, 23, 
    13, 15, 24, 
    16, 18, 25, 
    19, 21, 26, 
    27, 28, 29, 
    30, 31, 32, 
    33, 34, 35, 
    36, 37, 38, 
    39, 40, 41, 
    42, 43, 44, 
    45, 46, 42, 
    41, 47, 45, 
    42, 46, 43, 
    46, 48, 43, 
    45, 47, 46, 
    47, 49, 46, 
    46, 49, 48, 
    49, 38, 48, 
    41, 40, 47, 
    40, 50, 47, 
    47, 50, 49, 
    50, 51, 49, 
    49, 51, 38, 
    51, 36, 38, 
    43, 52, 44, 
    48, 53, 43, 
    38, 54, 48, 
    43, 53, 52, 
    53, 55, 52, 
    48, 54, 53, 
    54, 56, 53, 
    53, 56, 55, 
    56, 35, 55, 
    38, 37, 54, 
    37, 57, 54, 
    54, 57, 56, 
    57, 58, 56, 
    56, 58, 35, 
    58, 33, 35, 
    52, 59, 44, 
    55, 60, 52, 
    35, 61, 55, 
    52, 60, 59, 
    60, 62, 59, 
    55, 61, 60, 
    61, 63, 60, 
    60, 63, 62, 
    63, 32, 62, 
    35, 34, 61, 
    34, 64, 61, 
    61, 64, 63, 
    64, 65, 63, 
    63, 65, 32, 
    65, 30, 32, 
    59, 66, 44, 
    62, 67, 59, 
    32, 68, 62, 
    59, 67, 66, 
    67, 69, 66, 
    62, 68, 67, 
    68, 70, 67, 
    67, 70, 69, 
    70, 29, 69, 
    32, 31, 68, 
    31, 71, 68, 
    68, 71, 70, 
    71, 72, 70, 
    70, 72, 29, 
    72, 27, 29, 
    66, 42, 44, 
    69, 73, 66, 
    29, 74, 69, 
    66, 73, 42, 
    73, 45, 42, 
    69, 74, 73, 
    74, 75, 73, 
    73, 75, 45, 
    75, 41, 45, 
    29, 28, 74, 
    28, 76, 74, 
    74, 76, 75, 
    76, 77, 75, 
    75, 77, 41, 
    77, 39, 41, 
    78, 40, 39, 
    79, 80, 78, 
    26, 81, 79, 
    78, 80, 40, 
    80, 50, 40, 
    79, 81, 80, 
    81, 82, 80, 
    80, 82, 50, 
    82, 51, 50, 
    26, 21, 81, 
    21, 83, 81, 
    81, 83, 82, 
    83, 84, 82, 
    82, 84, 51, 
    84, 36, 51, 
    85, 37, 36, 
    86, 87, 85, 
    25, 88, 86, 
    85, 87, 37, 
    87, 57, 37, 
    86, 88, 87, 
    88, 89, 87, 
    87, 89, 57, 
    89, 58, 57, 
    25, 18, 88, 
    18, 90, 88, 
    88, 90, 89, 
    90, 91, 89, 
    89, 91, 58, 
    91, 33, 58, 
    92, 34, 33, 
    93, 94, 92, 
    24, 95, 93, 
    92, 94, 34, 
    94, 64, 34, 
    93, 95, 94, 
    95, 96, 94, 
    94, 96, 64, 
    96, 65, 64, 
    24, 15, 95, 
    15, 97, 95, 
    95, 97, 96, 
    97, 98, 96, 
    96, 98, 65, 
    98, 30, 65, 
    99, 31, 30, 
    100, 101, 99, 
    23, 102, 100, 
    99, 101, 31, 
    101, 71, 31, 
    100, 102, 101, 
    102, 103, 101, 
    101, 103, 71, 
    103, 72, 71, 
    23, 12, 102, 
    12, 104, 102, 
    102, 104, 103, 
    104, 105, 103, 
    103, 105, 72, 
    105, 27, 72, 
    106, 28, 27, 
    107, 108, 106, 
    22, 109, 107, 
    106, 108, 28, 
    108, 76, 28, 
    107, 109, 108, 
    109, 110, 108, 
    108, 110, 76, 
    110, 77, 76, 
    22, 9, 109, 
    9, 111, 109, 
    109, 111, 110, 
    111, 112, 110, 
    110, 112, 77, 
    112, 39, 77, 
    84, 85, 36, 
    83, 113, 84, 
    21, 114, 83, 
    84, 113, 85, 
    113, 86, 85, 
    83, 114, 113, 
    114, 115, 113, 
    113, 115, 86, 
    115, 25, 86, 
    21, 20, 114, 
    20, 116, 114, 
    114, 116, 115, 
    116, 117, 115, 
    115, 117, 25, 
    117, 16, 25, 
    91, 92, 33, 
    90, 118, 91, 
    18, 119, 90, 
    91, 118, 92, 
    118, 93, 92, 
    90, 119, 118, 
    119, 120, 118, 
    118, 120, 93, 
    120, 24, 93, 
    18, 17, 119, 
    17, 121, 119, 
    119, 121, 120, 
    121, 122, 120, 
    120, 122, 24, 
    122, 13, 24, 
    98, 99, 30, 
    97, 123, 98, 
    15, 124, 97, 
    98, 123, 99, 
    123, 100, 99, 
    97, 124, 123, 
    124, 125, 123, 
    123, 125, 100, 
    125, 23, 100, 
    15, 14, 124, 
    14, 126, 124, 
    124, 126, 125, 
    126, 127, 125, 
    125, 127, 23, 
    127, 10, 23, 
    105, 106, 27, 
    104, 128, 105, 
    12, 129, 104, 
    105, 128, 106, 
    128, 107, 106, 
    104, 129, 128, 
    129, 130, 128, 
    128, 130, 107, 
    130, 22, 107, 
    12, 11, 129, 
    11, 131, 129, 
    129, 131, 130, 
    131, 132, 130, 
    130, 132, 22, 
    132, 3, 22, 
    112, 78, 39, 
    111, 133, 112, 
    9, 134, 111, 
    112, 133, 78, 
    133, 79, 78, 
    111, 134, 133, 
    134, 135, 133, 
    133, 135, 79, 
    135, 26, 79, 
    9, 5, 134, 
    5, 136, 134, 
    134, 136, 135, 
    136, 137, 135, 
    135, 137, 26, 
    137, 19, 26, 
    138, 20, 19, 
    139, 140, 138, 
    8, 141, 139, 
    138, 140, 20, 
    140, 116, 20, 
    139, 141, 140, 
    141, 142, 140, 
    140, 142, 116, 
    142, 117, 116, 
    8, 7, 141, 
    7, 143, 141, 
    141, 143, 142, 
    143, 144, 142, 
    142, 144, 117, 
    144, 16, 117, 
    144, 17, 16, 
    143, 145, 144, 
    7, 146, 143, 
    144, 145, 17, 
    145, 121, 17, 
    143, 146, 145, 
    146, 147, 145, 
    145, 147, 121, 
    147, 122, 121, 
    7, 6, 146, 
    6, 148, 146, 
    146, 148, 147, 
    148, 149, 147, 
    147, 149, 122, 
    149, 13, 122, 
    149, 14, 13, 
    148, 150, 149, 
    6, 151, 148, 
    149, 150, 14, 
    150, 126, 14, 
    148, 151, 150, 
    151, 152, 150, 
    150, 152, 126, 
    152, 127, 126, 
    6, 2, 151, 
    2, 153, 151, 
    151, 153, 152, 
    153, 154, 152, 
    152, 154, 127, 
    154, 10, 127, 
    137, 138, 19, 
    136, 155, 137, 
    5, 156, 136, 
    137, 155, 138, 
    155, 139, 138, 
    136, 156, 155, 
    156, 157, 155, 
    155, 157, 139, 
    157, 8, 139, 
    5, 4, 156, 
    4, 158, 156, 
    156, 158, 157, 
    158, 1, 157, 
    157, 1, 8, 
    1, 0, 8, 
    154, 11, 10, 
    153, 159, 154, 
    2, 160, 153, 
    154, 159, 11, 
    159, 131, 11, 
    153, 160, 159, 
    160, 161, 159, 
    159, 161, 131, 
    161, 132, 131, 
    2, 1, 160, 
    1, 158, 160, 
    160, 158, 161, 
    158, 4, 161, 
    161, 4, 132, 
    4, 3, 132
};



#endif