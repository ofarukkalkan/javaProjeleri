package com.mycompany.moda;

/**
 *
 * @author omerfaruk
 */
public class moda {

    public static void main(String[] args) {
        // dim=5
        int dim = 5;
        // obj_no=2
        int obj_no = 2;
        // lb=0;
        double lb[][] = zerosMatrix(1, 1);
        // ub=1;
        double ub[][] = onesMatrix(1, 1);
        // if size(ub,2)==1
        if (getSizeOfDim(ub, 2) == 1) {
            // ub=ones(1,dim)*ub;
            ub = elementWiseProduct(onesMatrix(1, dim), ub[0][0]);
            // lb=ones(1,dim)*lb;
            lb = elementWiseProduct(onesMatrix(1, dim), lb[0][0]);
        }
        // moda algoritmasinin baslangic degerleri
        int max_iter = 100;
        int N = 100;
        int ArchiveMaxSize = 100;
        // Archive_X=zeros(100,dim);
        double Archive_X[][] = zerosMatrix(100, dim);
        // Archive_F=ones(100,obj_no)*inf;
        double Archive_F[][] = elementWiseProduct(onesMatrix(100, dim), Double.POSITIVE_INFINITY);
        //Archive_member_no = 0
        int Archive_member_no = 0;
        // r=(ub-lb)/2;
        double r[][] = elementWiseDivision(elementWiseSubstraction(ub, lb), 2);
        // V_max=(ub(1)-lb(1))/10;
        double V_max = (ub[0][1] - lb[0][1]) / 10.0;
        // Food_fitness = inf * ones(1, obj_no);
        double Food_fitness[] = elementWiseProduct(onesRowVector(obj_no), Double.POSITIVE_INFINITY);
        // Food_pos = zeros(dim, 1);
        double Food_pos[][] = zerosColumnVector(dim);
        // Enemy_fitness = -inf * ones(1, obj_no);
        double Enemy_fitness[] = elementWiseProduct(onesRowVector(obj_no), Double.NEGATIVE_INFINITY);
        // Enemy_pos = zeros(dim, 1);
        double Enemy_pos[][] = zerosColumnVector(dim);
        // X = initialization(N, dim, ub, lb);
        double X[][] = initialization(N, dim, ub[0], lb[0]);
        // fitness = zeros(N, 2);
        double fitness[][] = zerosMatrix(N, 2);
        // DeltaX = initialization(N, dim, ub, lb);
        double DeltaX[][] = initialization(N, dim, ub[0], lb[0]);
        //iter = 0;
        int iter = 0;
        // position_history = zeros(N, max_iter, dim);
        double position_history[][][] = zerosMatrix(N, max_iter, dim);
        double my_c, w, s, a, c, f, e;
        // max_iter = 100
        for (iter = 1; iter <= max_iter; iter++) {
            // r=(ub-lb)/4+((ub-lb)*(iter/max_iter)*2);
            r = elementWiseAddition(
                    elementWiseDivision(elementWiseSubstraction(ub, lb), 4),
                    elementWiseProduct(elementWiseSubstraction(ub, lb), (iter / max_iter) * 2)
            );
            // w=0.9-iter*((0.9-0.2)/max_iter);
            w = 0.9 - iter * ((0.9 - 0.2) / max_iter);
            // my_c=0.1-iter*((0.1-0)/(max_iter/2));
            my_c = 0.1 - iter * ((0.1 - 0) / (max_iter / 2));
            // if my_c<0
            if (iter < 0) {
                my_c = 0;
            }
            // if iter<(3*max_iter/4)
            if (iter < (3 * max_iter / 4)) {
                s = my_c;             // Seperation weight
                a = my_c;             // Alignment weight
                c = my_c;             // Cohesion weight
                f = 2 * Math.random();  // Food attraction weight
                e = my_c;             // Enemy distraction weight
            } else {
                s = my_c / iter;        // Seperation weight
                a = my_c / iter;        // Alignment weight
                c = my_c / iter;        // Coheikomuksdfssion weight
                f = 2 * Math.random();  // Food attraction weight
                e = my_c / iter;        // Enemy distraction weight
            }
            // Particles_F=zeros(1,2);
            double Particles_F[][] = new double[N][obj_no];
            // for i=1:N Calculate all the objective values first
            for (int i = 0; i < N; i++) {
                double tmp[] = ZDT1(transpose(getJthColumn(X, i)));
                // Particles_F(i,:) =ObjectiveFunction(X(:,i)');
                for (int j = 0; j < obj_no; j++) {
                    Particles_F[i][j] = tmp[i];
                    Particles_F[i][j] = tmp[i];
                }
                // if dominates(Particles_F(i,:),Food_fitness)
                if (dominates(Particles_F[i], Food_fitness)) {
                    // Food_fitness=Particles_F(i,:);
                    Food_fitness = getIthRow(Particles_F, i);
                    // Food_pos=X(:,i);
                    Food_pos = getJthColumn(X, i);
                }
                // if dominates(Enemy_fitness,Particles_F(i,:))
                if (dominates(Enemy_fitness, Particles_F[i])) {
                    // if all(X(:,i)<ub') && all( X(:,i)>lb')
                    if ( // all(X(:,i)<ub')
                            isAllZero(isLessThen(getJthColumn(X, i), 0, transpose(ub), 0))
                            && // all( X(:,i)>lb')
                            isAllZero(isGreaterThen(getJthColumn(X, i), 0, transpose(lb), 0))) {

                        // Enemy_fitness=Particles_F(i,:);
                        Enemy_fitness = getIthRow(Particles_F, i);
                        // Enemy_pos=X(:,i);
                        Enemy_pos = getJthColumn(X, i);
                    }
                }
            }
            // [Archive_X, Archive_F, Archive_member_no]=UpdateArchive(Archive_X, Archive_F, X, Particles_F, Archive_member_no);
            Object tmp[] = new Object[3];
            tmp = UpdateArchive(Archive_X, Archive_F, Particles_F, Particles_F, Archive_member_no);
            Archive_X = (double[][])tmp[0];
            Archive_F = (double[][])tmp[1];
            Archive_member_no = (int)tmp[2];
            
            // if Archive_member_no>ArchiveMaxSize
            if(Archive_member_no > ArchiveMaxSize){
                // RankingProcess yazilacakti en son !!
            }
            
        }
    }

    public static Object[] UpdateArchive(double Archive_X[][], double Archive_F[][], double Particles_X[][], double Particles_F[][], int Archive_member_no) {
        Object results[] = new Object[3];
        // Archive_X_temp=[Archive_X ; Particles_X'];
        double Archive_X_temp[][] = matrixFromMatrices(Archive_X, transpose(Archive_X));
        // Archive_F_temp=[Archive_F ; Particles_F];
        double Archive_F_temp[][] = matrixFromMatrices(Archive_F, Particles_F);
        // o=zeros(1,size(Archive_F_temp,1));
        double o[] = zerosRowVector(getSizeOfDim(Archive_F_temp, 1));
        // for i=1:size(Archive_F_temp,1)
        for (int i = 0; i < getSizeOfDim(Archive_F_temp, 1); i++) {
            // o(i)=0;
            o[i] = 0;
            // for j=1:i-1
            for (int j = 0; j < i - 1; j++) {
                // if any(Archive_F_temp(i,:) ~= Archive_F_temp(j,:))
                if (isAnyZero(isNotEqual(getIthRow(Archive_F_temp, i), getIthRow(Archive_F_temp, j)))) {
                    // if dominates(Archive_F_temp(i,:),Archive_F_temp(j,:))
                    if (dominates(getIthRow(Archive_F_temp, i), getIthRow(Archive_F_temp, j))) {
                        // o(j)=1;
                        o[j] = 1;
                        // elseif dominates(Archive_F_temp(j,:),Archive_F_temp(i,:))    
                    } else if (dominates(getIthRow(Archive_F_temp, j), getIthRow(Archive_F_temp, i))) {
                        // o(i)=1;
                        o[i] = 1;
                        // break;
                        break;
                    }
                } else {
                    // o(j)=1;
                    // o(i)=1;
                    o[j] = 1;
                    o[i] = 1;
                }
            }
        }
        // Archive_member_no=0;
        // index=0;

        Archive_member_no = 0;
        int index = 0;
        double Archive_X_updated[][] = new double[Archive_X.length][Archive_X[0].length];
        double Archive_F_updated[][] = new double[Archive_F.length][Archive_F[0].length];

        // for i=1:size(Archive_X_temp,1)
        for (int i = 0; i < getSizeOfDim(Archive_X_temp, 1); i++) {
            // if o(i)==0
            if (o[i] == 0) {

                // Archive_X_updated(Archive_member_no,:)=Archive_X_temp(i,:);
                for (int j = 0; j < Archive_X_updated.length; j++) {
                    Archive_X_updated[Archive_member_no] = getIthRow(Archive_X_temp, i);
                }
                // Archive_F_updated(Archive_member_no,:)=Archive_F_temp(i,:);
                for (int j = 0; j < Archive_F_updated.length; j++) {
                    Archive_F_updated[Archive_member_no] = getIthRow(Archive_F_temp, i);
                }
                // Archive_member_no=Archive_member_no+1;
                Archive_member_no++;
            } else {
                index++;
            }
        }
        results[0] = Archive_X_updated;
        results[1] = Archive_F_updated;
        results[2] = Archive_member_no;

        return results;
    }

    public static boolean dominates(double x[], double y[]) {
        return isAllZero(isLessThenAndEqual(x, y)) && isAnyZero(isLessThen(x, y));

    }

    public static double[] ZDT1(double arr[][]) {

        return arr[0];
    }

    public static double[][] initialization(int SearchAgents_no, int dim, double ub[], double lb[]) {
        // Boundary_no= size(ub,2);
        int Boundary_no = ub.length;

        double ub_new[] = new double[dim];
        double lb_new[] = new double[dim];

        if (Boundary_no == 1) {
            //ub_new=ones(1,dim)*ub;
            for (int i = 0; i < dim; i++) {
                ub_new[i] = ub[0];
            }
            //lb_new=ones(1,dim)*lb;
            for (int i = 0; i < dim; i++) {
                lb_new[i] = lb[0];
            }
        } else {
            // ub_new=ub;       
            ub_new = ub;
            // lb_new=lb;
            lb_new = lb;
        }

        double X[][] = new double[SearchAgents_no][dim];

        /*
        for i=1:dim
        ub_i=ub_new(i);
        lb_i=lb_new(i);
        X(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
        end
         */
        for (int i = 0; i < dim; i++) {
            double ub_i = ub_new[i];
            double lb_i = lb_new[i];
            for (int j = 0; j < SearchAgents_no; j++) {
                X[j][i] = Math.random() * (ub_i - lb_i) + lb_i;
            }

        }

        // X=X';
        return transpose(X);
    }

    public static double[][] transpose(double[][] m) {
        double[][] temp = new double[m[0].length][m.length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                temp[j][i] = m[i][j];
            }
        }
        return temp;
    }

    public static double[] tanspose(double m[][]) {
        double arr[] = new double[m[0].length];
        for (int i = 0; i < m[0].length; i++) {
            arr[i] = m[i][0];
        }
        return arr;
    }

    public static double[][] tanspose(double vec[]) {
        double m[][] = new double[vec.length][1];
        for (int i = 0; i < vec.length; i++) {
            m[i][0] = vec[i];
        }
        return m;
    }

    public static double[] getIthRow(double matrix[][], int I) {
        double arr[] = new double[matrix[0].length];
        for (int i = 0; i < matrix[0].length; i++) {
            arr[i] = matrix[I][i];
        }
        return arr;
    }

    public static double[][] getJthColumn(double matrix[][], int J) {
        double arr[][] = new double[matrix.length][1];
        for (int i = 0; i < matrix.length; i++) {
            arr[i][0] = matrix[i][J];
        }
        return arr;
    }

    public static double[][] zerosMatrix(int r, int c) {
        double m[][] = new double[r][c];
        return m;
    }

    public static double[][] onesMatrix(int r, int c) {
        double m[][] = new double[r][c];
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                m[i][j] = 1;
            }
        }
        return m;
    }

    public static double[][][] zerosMatrix(int r, int c, int p) {
        double m[][][] = new double[r][c][p];
        return m;
    }

    public static double[][][] onesMatrix(int r, int c, int p) {
        double m[][][] = new double[r][c][p];
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                for (int k = 0; k < p; k++) {
                    m[i][j][k] = 1;
                }
            }
        }
        return m;
    }

    public static double[] zerosRowVector(int r) {
        double m[] = new double[r];
        return m;
    }

    public static double[] onesRowVector(int r) {
        double m[] = new double[r];
        for (int i = 0; i < r; i++) {
            m[i] = 1;
        }
        return m;
    }

    public static double[][] zerosColumnVector(int c) {
        double m[][] = new double[c][1];
        return m;
    }

    public static double[][] onesColumnVector(int c) {
        double m[][] = new double[c][1];
        for (int i = 0; i < c; i++) {
            m[i][0] = 1;
        }
        return m;
    }

    public static double[][] matrixFromMatrices(double m1[][], double m2[][]) {
        double m3[][] = new double[m1.length + m2.length][m1[0].length];
        for (int i = 0; i < m1.length; i++) {
            for (int j = 0; j < m1[0].length; j++) {
                m3[i][j] = m1[i][j];
            }
        }
        for (int i = m1.length; i < m2.length; i++) {
            for (int j = 0; j < m2[0].length; j++) {
                m3[i][j] = m2[i][j];
            }
        }
        return m3;
    }

    public static double[][] elementWiseDivision(double m[][], double scalar) {
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                m[i][j] = m[i][j] / scalar;
            }
        }
        return m;
    }

    public static double[][] elementWiseProduct(double m[][], double scalar) {
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                m[i][j] = m[i][j] * scalar;
            }
        }
        return m;
    }

    public static double[][] elementWiseAddition(double m[][], double scalar) {
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                m[i][j] = m[i][j] + scalar;
            }
        }
        return m;
    }

    public static double[][] elementWiseSubstraction(double m[][], double scalar) {
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                m[i][j] = m[i][j] - scalar;
            }
        }
        return m;
    }

    public static double[][] elementWiseSubstraction(double m[][], double n[][]) {
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                m[i][j] = m[i][j] - n[i][j];
            }
        }
        return m;
    }

    public static double[][] elementWiseAddition(double m[][], double n[][]) {
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                m[i][j] = m[i][j] + n[i][j];
            }
        }
        return m;
    }

    public static double[] elementWiseProduct(double vec[], double scalar) {
        for (int i = 0; i < vec.length; i++) {
            vec[i] *= scalar;
        }
        return vec;
    }

    public static int getSizeOfDim(double m[], int dim) {
        int size = 0;
        switch (dim) {
            case 1:
                size = m.length;
                break;
        }
        return size;
    }

    public static int getSizeOfDim(double m[][], int dim) {
        int size = 0;
        switch (dim) {
            case 1:
                size = m.length;
                break;
            case 2:
                size = m[0].length;
                break;
        }
        return size;
    }

    public static int getSizeOfDim(double m[][][], int dim) {
        int size = 0;
        switch (dim) {
            case 1:
                size = m.length;
                break;
            case 2:
                size = m[0].length;
                break;
            case 3:
                size = m[0][0].length;
                break;
        }
        return size;
    }

    public static boolean[] isAllZero(double m[][]) {
        boolean result[] = new boolean[m.length];
        for (int i = 0; i < m[0].length; i++) {
            result[i] = isAllZero(m, i);
        }
        return result;
    }

    public static boolean isAllZero(double m[]) {
        for (int i = 0; i < m.length; i++) {
            if (m[i] != 0) {
                return false;
            }
        }
        return true;
    }

    public static boolean isAllZero(double m[][], int c) {
        for (int i = 0; i < m.length; i++) {
            if (m[i][c] != 0) {
                return false;
            }
        }
        return true;
    }

    public static boolean[] isAnyZero(double m[][]) {
        boolean result[] = new boolean[m.length];
        for (int i = 0; i < m[0].length; i++) {
            result[i] = isAnyZero(m, i);
        }
        return result;
    }

    public static boolean isAnyZero(double m[]) {
        for (int i = 0; i < m.length; i++) {
            if (m[i] != 0) {
                return false;
            }
        }
        return true;
    }

    public static boolean isAnyZero(double m[][], int c) {
        for (int i = 0; i < m.length; i++) {
            if (m[i][c] != 0) {
                return false;
            }
        }
        return true;
    }

    public static double[] isLessThenAndEqual(double m1[], double m2[]) {
        double res[] = new double[m1.length];
        for (int i = 0; i < m1.length; i++) {
            if (m1[i] <= m2[i]) {
                res[i] = 1;
            } else {
                res[i] = 0;
            }
        }
        return res;
    }

    public static double[] isLessThen(double m1[][], int c1, double m2[][], int c2) {
        double res[] = new double[m1.length];
        for (int i = 0; i < m1.length; i++) {
            if (m1[i][c1] < m2[i][c2]) {
                res[i] = 1;
            } else {
                res[i] = 0;
            }
        }
        return res;
    }

    public static double[] isLessThen(double m1[], double m2[]) {
        double res[] = new double[m1.length];
        for (int i = 0; i < m1.length; i++) {
            if (m1[i] < m2[i]) {
                res[i] = 1;
            } else {
                res[i] = 0;
            }
        }
        return res;
    }

    public static double[] isGreaterThen(double m1[][], int c1, double m2[][], int c2) {
        double res[] = new double[m1.length];
        for (int i = 0; i < m1.length; i++) {
            if (m1[i][c1] > m2[i][c2]) {
                res[i] = 1;
            } else {
                res[i] = 0;
            }
        }
        return res;
    }

    public static double[] isNotEqual(double m1[], double m2[]) {
        double res[] = new double[m1.length];
        for (int i = 0; i < m1.length; i++) {
            if (m1[i] != m2[i]) {
                res[i] = 1;
            } else {
                res[i] = 0;
            }
        }
        return res;
    }

}
