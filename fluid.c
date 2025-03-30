#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <SDL2/SDL.h>

#define IX(x, y, N) ((x) + (y) * (N))

typedef struct FluidCube
{
    int size;
    float dt;
    float diff;
    float visc;

    float *s;
    float *density;

    float *Vx;
    float *Vy;

    float *Vx0;
    float *Vy0;
} FluidCube;

FluidCube *FluidCubeCreate(int size, float diffusion, float viscosity, float dt)
{
    FluidCube *cube = (FluidCube *)malloc(sizeof(*cube));

    int N = size;

    cube->size = size;
    cube->dt = dt;
    cube->diff = diffusion;
    cube->visc = viscosity;

    cube->s = calloc(N * N, sizeof(float));
    cube->density = calloc(N * N, sizeof(float));

    cube->Vx = calloc(N * N, sizeof(float));
    cube->Vy = calloc(N * N, sizeof(float));

    cube->Vx0 = calloc(N * N, sizeof(float));
    cube->Vy0 = calloc(N * N, sizeof(float));

    return cube;
}

void FluidCubeFree(FluidCube *cube)
{
    free(cube->s);
    free(cube->density);
    free(cube->Vx);
    free(cube->Vy);
    free(cube->Vx0);
    free(cube->Vy0);
    free(cube);
}

static void set_bnd(int b, float *x, int N)
{
    for (int i = 1; i < N - 1; i++)
    {
        x[IX(i, 0, N)] = (b == 2) ? -x[IX(i, 1, N)] : x[IX(i, 1, N)];
        x[IX(i, N - 1, N)] = (b == 2) ? -x[IX(i, N - 2, N)] : x[IX(i, N - 2, N)];
    }
    for (int j = 1; j < N - 1; j++)
    {
        x[IX(0, j, N)] = (b == 1) ? -x[IX(1, j, N)] : x[IX(1, j, N)];
        x[IX(N - 1, j, N)] = (b == 1) ? -x[IX(N - 2, j, N)] : x[IX(N - 2, j, N)];
    }

    x[IX(0, 0, N)] = 0.5f * (x[IX(1, 0, N)] + x[IX(0, 1, N)]);
    x[IX(0, N - 1, N)] = 0.5f * (x[IX(1, N - 1, N)] + x[IX(0, N - 2, N)]);
    x[IX(N - 1, 0, N)] = 0.5f * (x[IX(N - 2, 0, N)] + x[IX(N - 1, 1, N)]);
    x[IX(N - 1, N - 1, N)] = 0.5f * (x[IX(N - 2, N - 1, N)] + x[IX(N - 1, N - 2, N)]);
}

static void lin_solve(int b, float *x, float *x0, float a, float c, int iter, int N)
{
    float cRecip = 1.0f / c;

    for (int k = 0; k < iter; k++)
    {
        for (int j = 1; j < N - 1; j++)
        {
            for (int i = 1; i < N - 1; i++)
            {
                x[IX(i, j, N)] = (x0[IX(i, j, N)] + a * (x[IX(i + 1, j, N)] + x[IX(i - 1, j, N)] + x[IX(i, j + 1, N)] + x[IX(i, j - 1, N)])) * cRecip;
            }
        }

        set_bnd(b, x, N);
    }
}

static void diffuse(int b, float *x, float *x0, float diff, float dt, int iter, int N)
{
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 4 * a, iter, N);
}

static void advect(int b, float *d, float *d0, float *velocX, float *velocY, float dt, int N)
{
    float i0, i1, j0, j1;
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    float s0, s1, t0, t1, tmp1, tmp2, x, y;
    float Nfloat = (float)N;

    for (int j = 1; j < N - 1; j++)
    {
        for (int i = 1; i < N - 1; i++)
        {
            tmp1 = dtx * velocX[IX(i, j, N)];
            tmp2 = dty * velocY[IX(i, j, N)];

            x = i - tmp1;
            y = j - tmp2;

            if (x < 0.5f) x = 0.5f;
            if (x > Nfloat - 0.5f) x = Nfloat - 0.5f;

            i0 = floorf(x);
            i1 = i0 + 1;

            if (y < 0.5f) y = 0.5f;
            if (y > Nfloat - 0.5f) y = Nfloat - 0.5f;

            j0 = floorf(y);
            j1 = j0 + 1;
            s1 = x - i0;
            s0 = 1.0f - s1;
            t1 = y - j0;
            t0 = 1.0f - t1;

            int i0i = (int)i0, i1i = (int)i1;
            int j0i = (int)j0, j1i = (int)j1;

            d[IX(i, j, N)] = s0 * (t0 * d0[IX(i0i, j0i, N)] + t1 * d0[IX(i0i, j1i, N)]) + s1 * (t0 * d0[IX(i1i, j0i, N)] + t1 * d0[IX(i1i, j1i, N)]);
        }
    }

    set_bnd(b, d, N);
}

static void project(float *velocX, float *velocY, float *p, float *div, int iter, int N)
{
    for (int j = 1; j < N - 1; j++)
    {
        for (int i = 1; i < N - 1; i++)
        {
            div[IX(i, j, N)] = -0.5f * (velocX[IX(i + 1, j, N)] - velocX[IX(i - 1, j, N)] + velocY[IX(i, j + 1, N)] - velocY[IX(i, j - 1, N)]) / N;
            p[IX(i, j, N)] = 0;
        }
    }

    set_bnd(0, div, N);
    set_bnd(0, p, N);
    lin_solve(0, p, div, 1, 4, iter, N);

    for (int j = 1; j < N - 1; j++)
    {
        for (int i = 1; i < N - 1; i++)
        {
            velocX[IX(i, j, N)] -= 0.5f * (p[IX(i + 1, j, N)] - p[IX(i - 1, j, N)]) * N;
            velocY[IX(i, j, N)] -= 0.5f * (p[IX(i, j + 1, N)] - p[IX(i, j - 1, N)]) * N;
        }
    }

    set_bnd(1, velocX, N);
    set_bnd(2, velocY, N);
}

void FluidCubeStep(FluidCube *cube)
{
    int N = cube->size;
    float visc = cube->visc;
    float diff = cube->diff;
    float dt = cube->dt;
    float *Vx = cube->Vx;
    float *Vy = cube->Vy;
    float *Vx0 = cube->Vx0;
    float *Vy0 = cube->Vy0;
    float *s = cube->s;
    float *density = cube->density;

    diffuse(1, Vx0, Vx, visc, dt, 4, N);
    diffuse(2, Vy0, Vy, visc, dt, 4, N);

    float *p = calloc(N * N, sizeof(float));
    float *div = calloc(N * N, sizeof(float));

    project(Vx0, Vy0, p, div, 4, N);

    advect(1, Vx, Vx0, Vx0, Vy0, dt, N);
    advect(2, Vy, Vy0, Vx0, Vy0, dt, N);

    project(Vx, Vy, p, div, 4, N);

    diffuse(0, s, density, diff, dt, 4, N);
    advect(0, density, s, Vx, Vy, dt, N);

    free(p);
    free(div);
}

void FluidCubeAddDensity(FluidCube *cube, int x, int y, float amount)
{
    int N = cube->size;
    cube->density[IX(x, y, N)] += amount;
}

void FluidCubeAddVelocity(FluidCube *cube, int x, int y, float amountX, float amountY)
{
    int N = cube->size;
    int index = IX(x, y, N);
    cube->Vx[index] += amountX;
    cube->Vy[index] += amountY;
}

int main(int argc, char *argv[])
{
    const int gridSize = 64;
    const float dt = 0.1f;
    const float diffusion = 0.0001f;
    const float viscosity = 0.0001f;
    const int windowSize = 512;

    FluidCube *cube = FluidCubeCreate(gridSize, diffusion, viscosity, dt);

    SDL_Window *window = SDL_CreateWindow("Fluid Simulation", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, windowSize, windowSize, SDL_WINDOW_SHOWN);

    SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    int running = 1;

    SDL_Event event;

    int mouseDown = 0;
    int prevMouseX = 0, prevMouseY = 0;

    while (running)
    {
        while (SDL_PollEvent(&event))
        {
            if (event.type == SDL_QUIT) running = 0;
            if (event.type == SDL_MOUSEBUTTONDOWN)
            {
                if (event.button.button == SDL_BUTTON_LEFT)
                {
                    mouseDown = 1;
                    prevMouseX = event.button.x;
                    prevMouseY = event.button.y;
                }
            }
            if (event.type == SDL_MOUSEBUTTONUP)
            {
                if (event.button.button == SDL_BUTTON_LEFT)
                {
                    mouseDown = 0;
                }
            }
        }

        if (mouseDown)
        {
            int mouseX, mouseY;
            SDL_GetMouseState(&mouseX, &mouseY);
            
            int gridX = mouseX * gridSize / windowSize;
            int gridY = mouseY * gridSize / windowSize;

            
            FluidCubeAddDensity(cube, gridX, gridY, 100.0f);
            
            float forceX = (mouseX - prevMouseX) * 0.1f;
            float forceY = (mouseY - prevMouseY) * 0.1f;

            FluidCubeAddVelocity(cube, gridX, gridY, forceX, forceY);

            prevMouseX = mouseX;
            prevMouseY = mouseY;
        }

        
        FluidCubeStep(cube);

        
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        for (int j = 0; j < gridSize; j++)
        {
            for (int i = 0; i < gridSize; i++)
            {
                float d = cube->density[IX(i, j, gridSize)];
                
                int intensity = (int)(d * 10.0f);
                if (intensity > 255) intensity = 255;
                if (intensity < 0) intensity = 0;

                Uint8 r = (Uint8)intensity;
                Uint8 g = 0;
                Uint8 b = (Uint8)(255 - intensity);

                SDL_SetRenderDrawColor(renderer, r, g, b, 255);

                SDL_Rect rect = {i * (windowSize / gridSize), j * (windowSize / gridSize), windowSize / gridSize, windowSize / gridSize};

                SDL_RenderFillRect(renderer, &rect);
            }
        }

        SDL_RenderPresent(renderer);
        SDL_Delay(16); // frame time
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    FluidCubeFree(cube);
    return 0;
}
