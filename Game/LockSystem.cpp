#include "LockSystem.h"

LockSystem::LockSystem() {
    wallLock.fill(false);
    wallAngle.fill(0);
}
void LockSystem::lockWalls(int indexLayer) {
    switch (indexLayer) {
    case FRONT:
    case BACK:
        wallLock[LEFT] = wallLock[RIGHT] = wallLock[UP] = wallLock[DOWN] = true;
        break;
    case LEFT:
    case RIGHT:
        wallLock[FRONT] = wallLock[BACK] = wallLock[UP] = wallLock[DOWN] = true;
        break;
    case UP:
    case DOWN:
        wallLock[FRONT] = wallLock[BACK] = wallLock[LEFT] = wallLock[RIGHT] = true;
        break;
    }
}

void LockSystem::unlockWalls(int indexLayer) {
    switch (indexLayer) {
    case FRONT:
        if (wallAngle[BACK] % 90 == 0) wallLock[LEFT] = wallLock[RIGHT] = wallLock[UP] = wallLock[DOWN] = false;
        break;
    case BACK:
        if (wallAngle[FRONT] % 90 == 0) wallLock[LEFT] = wallLock[RIGHT] = wallLock[UP] = wallLock[DOWN] = false;
        break;
    case LEFT:
        if (wallAngle[RIGHT] % 90 == 0) wallLock[FRONT] = wallLock[BACK] = wallLock[UP] = wallLock[DOWN] = false;
        break;
    case RIGHT:
        if (wallAngle[LEFT] % 90 == 0) wallLock[FRONT] = wallLock[BACK] = wallLock[UP] = wallLock[DOWN] = false;
        break;
    case UP:
        if (wallAngle[DOWN] % 90 == 0) wallLock[FRONT] = wallLock[BACK] = wallLock[LEFT] = wallLock[RIGHT] = false;
        break;
    case DOWN:
        if (wallAngle[UP] % 90 == 0) wallLock[FRONT] = wallLock[BACK] = wallLock[LEFT] = wallLock[RIGHT] = false;
        break;
    }
}