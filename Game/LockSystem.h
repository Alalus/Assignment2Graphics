// LockSystem.h
#pragma once
#include <array>

// Define an enumeration for the cube faces
enum CubeFace {
    FRONT,
    BACK,
    LEFT,
    RIGHT,
    UP,
    DOWN
};

class LockSystem {
private:
    std::array<bool, 6> wallLock; // Tracks whether each face is locked
    std::array<int, 6> wallAngle; // Tracks the rotation angle of each face

public:
    LockSystem();
    void lockWalls(int indexLayer);
    void unlockWalls(int indexLayer);
    bool isWallLocked(int indexLayer) const {
        return wallLock[indexLayer];
    }
    int getWallAngel(int indexLayer) const {
        return wallAngle[indexLayer];
    }
    void setWallAngel(int angle, int indexLayer) {
        wallAngle[indexLayer] = angle;
    }
    void modulusAngel(int indexLayer) {
        wallAngle[indexLayer] = wallAngle[indexLayer] % 360;
    }

};