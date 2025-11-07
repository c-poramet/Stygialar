// Light Out Puzzle Solver Module
const LightOut = {
    gridWidth: 5,
    gridHeight: 5,
    maxStates: 2,
    cellStates: [],
    stateConfig: null, // null = uniform states, array = custom states per cell

    reset() {
        this.cellStates = new Array(this.gridWidth * this.gridHeight).fill(0);
        this.stateConfig = null;
    },

    initialize(width, height, states) {
        this.gridWidth = width;
        this.gridHeight = height;
        this.maxStates = states;
        this.reset();
    },

    setState(states) {
        if (states.length !== this.gridHeight * this.gridWidth) {
            throw new Error(`Expected ${this.gridHeight * this.gridWidth} states, got ${states.length}`);
        }
        
        // Validate states
        for (let i = 0; i < states.length; i++) {
            const maxState = this.stateConfig ? this.stateConfig[i] : this.maxStates;
            if (states[i] < 0 || states[i] >= maxState) {
                throw new Error(`Invalid state ${states[i]} at position ${i}. Must be 0-${maxState - 1}`);
            }
        }
        
        this.cellStates = [...states];
    },

    setCustomStateConfig(config) {
        if (config.length !== this.gridHeight * this.gridWidth) {
            throw new Error(`Expected ${this.gridHeight * this.gridWidth} state configs, got ${config.length}`);
        }
        this.stateConfig = [...config];
    },

    getState() {
        return [...this.cellStates];
    },

    getStateAsGrid() {
        const grid = [];
        for (let i = 0; i < this.gridHeight; i++) {
            const row = this.cellStates.slice(
                i * this.gridWidth,
                (i + 1) * this.gridWidth
            );
            grid.push(row);
        }
        return grid;
    },

    // Modular inverse using extended Euclidean algorithm
    modInverse(a, m) {
        a = ((a % m) + m) % m;
        for (let x = 1; x < m; x++) {
            if ((a * x) % m === 1) {
                return x;
            }
        }
        return null;
    },

    // Toggle a cell and its neighbors
    toggleCell(row, col, times = 1) {
        const index = row * this.gridWidth + col;
        const maxState = this.stateConfig ? this.stateConfig[index] : this.maxStates;
        
        // Toggle center
        this.cellStates[index] = (this.cellStates[index] + times) % maxState;
        
        // Toggle neighbors
        const neighbors = [
            [row - 1, col],     // up
            [row + 1, col],     // down
            [row, col - 1],     // left
            [row, col + 1]      // right
        ];
        
        for (const [r, c] of neighbors) {
            if (r >= 0 && r < this.gridHeight && c >= 0 && c < this.gridWidth) {
                const nIndex = r * this.gridWidth + c;
                const nMaxState = this.stateConfig ? this.stateConfig[nIndex] : this.maxStates;
                this.cellStates[nIndex] = (this.cellStates[nIndex] + times) % nMaxState;
            }
        }
    },

    // Solve using Gaussian elimination over GF(n)
    solve() {
        const width = this.gridWidth;
        const height = this.gridHeight;
        const totalCells = width * height;
        
        // Build augmented matrix
        const matrix = [];
        
        for (let i = 0; i < totalCells; i++) {
            const row = new Array(totalCells + 1).fill(0);
            const cellRow = Math.floor(i / width);
            const cellCol = i % width;
            const maxState = this.stateConfig ? this.stateConfig[i] : this.maxStates;
            
            // Cell affects itself
            row[i] = 1;
            
            // Cell affects neighbors
            if (cellRow > 0) {
                row[(cellRow - 1) * width + cellCol] = 1; // up
            }
            if (cellRow < height - 1) {
                row[(cellRow + 1) * width + cellCol] = 1; // down
            }
            if (cellCol > 0) {
                row[cellRow * width + (cellCol - 1)] = 1; // left
            }
            if (cellCol < width - 1) {
                row[cellRow * width + (cellCol + 1)] = 1; // right
            }
            
            // Target state (want to reach 0 from current state)
            row[totalCells] = (maxState - this.cellStates[i]) % maxState;
            
            matrix.push(row);
        }
        
        // Gaussian elimination in GF(maxStates)
        const mod = this.maxStates;
        let pivot = 0;
        
        for (let col = 0; col < totalCells && pivot < totalCells; col++) {
            // Find pivot
            let pivotRow = -1;
            for (let row = pivot; row < totalCells; row++) {
                if (matrix[row][col] % mod !== 0) {
                    pivotRow = row;
                    break;
                }
            }
            
            if (pivotRow === -1) continue;
            
            // Swap rows
            [matrix[pivot], matrix[pivotRow]] = [matrix[pivotRow], matrix[pivot]];
            
            // Make pivot element 1
            const pivotVal = matrix[pivot][col];
            const pivotInv = this.modInverse(pivotVal, mod);
            if (pivotInv === null) continue;
            
            for (let j = 0; j <= totalCells; j++) {
                matrix[pivot][j] = (matrix[pivot][j] * pivotInv) % mod;
            }
            
            // Eliminate column
            for (let row = 0; row < totalCells; row++) {
                if (row === pivot) continue;
                const factor = matrix[row][col];
                for (let j = 0; j <= totalCells; j++) {
                    matrix[row][j] = (matrix[row][j] - factor * matrix[pivot][j] % mod + mod) % mod;
                }
            }
            
            pivot++;
        }
        
        // Check for inconsistency
        for (let i = pivot; i < totalCells; i++) {
            let allZero = true;
            for (let j = 0; j < totalCells; j++) {
                if (matrix[i][j] % mod !== 0) {
                    allZero = false;
                    break;
                }
            }
            if (allZero && matrix[i][totalCells] % mod !== 0) {
                return { solvable: false, solution: null, message: 'No solution exists' };
            }
        }
        
        // Extract solution
        const solution = [];
        const solutionGrid = new Array(totalCells).fill(0);
        
        for (let i = 0; i < totalCells; i++) {
            // Find leading 1
            let leadingCol = -1;
            for (let j = 0; j < totalCells; j++) {
                if (matrix[i][j] % mod !== 0) {
                    leadingCol = j;
                    break;
                }
            }
            
            if (leadingCol !== -1) {
                const clicks = matrix[i][totalCells] % mod;
                if (clicks !== 0) {
                    solutionGrid[leadingCol] = clicks;
                    for (let k = 0; k < clicks; k++) {
                        solution.push(leadingCol);
                    }
                }
            }
        }
        
        return {
            solvable: true,
            solution,
            solutionGrid,
            clicks: solution.length,
            uniqueClicks: solutionGrid.filter(x => x > 0).length
        };
    },

    // Simulate solution with step-by-step states
    simulateSolution(solution) {
        const steps = [];
        const state = [...this.cellStates];
        
        // Record initial state
        steps.push({
            step: 0,
            action: 'Initial State',
            state: [...state]
        });
        
        // Apply each click
        const clickCounts = {};
        for (const cellIndex of solution) {
            clickCounts[cellIndex] = (clickCounts[cellIndex] || 0) + 1;
        }
        
        for (const [cellIndex, times] of Object.entries(clickCounts)) {
            const idx = parseInt(cellIndex);
            const row = Math.floor(idx / this.gridWidth);
            const col = idx % this.gridWidth;
            
            // Save current state
            const tempState = this.cellStates;
            this.cellStates = [...state];
            
            // Perform toggle
            this.toggleCell(row, col, times);
            
            // Record new state
            state.splice(0, state.length, ...this.cellStates);
            steps.push({
                step: steps.length,
                action: `Click (${row + 1}, ${col + 1}) ${times > 1 ? times + ' times' : ''}`,
                cellIndex: idx,
                row: row + 1,
                col: col + 1,
                times,
                state: [...state]
            });
            
            // Restore state
            this.cellStates = tempState;
        }
        
        return steps;
    }
};
