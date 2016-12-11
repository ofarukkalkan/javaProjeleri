package com.mycompany.moda;

/**
 *
 * @author omerfaruk
 */
public class moda {

    public static void main(String[] args) {
        moda();
    }

    public static void moda() {
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
        double Archive_F[][] = elementWiseProduct(onesMatrix(100, obj_no), Double.POSITIVE_INFINITY);
        //Archive_member_no = 0
        int Archive_member_no = 0;
        // r=(ub-lb)/2;
        double r[][] = elementWiseDivision(elementWiseSubstraction(ub, lb), 2.0);
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
                    elementWiseProduct(elementWiseSubstraction(ub, lb), (iter * 1.0f / max_iter) * 2)
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
            if (iter < (3.0 * max_iter / 4.0)) {
                s = my_c;               // Seperation weight
                a = my_c;               // Alignment weight
                c = my_c;               // Cohesion weight
                f = 2.0 * Math.random();// Food attraction weight
                e = my_c;               // Enemy distraction weight
            } else {
                s = my_c / iter;        // Seperation weight
                a = my_c / iter;        // Alignment weight
                c = my_c / iter;        // Coheikomuksdfssion weight
                f = 2.0 * Math.random();// Food attraction weight
                e = my_c / iter;        // Enemy distraction weight
            }
            // Particles_F=zeros(1,2);
            double Particles_F[][] = new double[N][obj_no];
            // for i=1:N Calculate all the objective values first
            for (int i = 0; i < N; i++) {
                // Particles_F(i,:) =ObjectiveFunction(X(:,i)');
                setRow(Particles_F, i, ZDT1(transposeColumn(getJthColumn(X, i))));
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
                            isAllNonZero(isLessThen(getJthColumn(X, i), 0, transposeMatrix(ub), 0))
                            && // all( X(:,i)>lb')
                            isAllNonZero(isGreaterThen(getJthColumn(X, i), 0, transposeMatrix(lb), 0))) {

                        // Enemy_fitness=Particles_F(i,:);
                        Enemy_fitness = getIthRow(Particles_F, i);
                        // Enemy_pos=X(:,i);
                        Enemy_pos = getJthColumn(X, i);
                    }
                }
            }
            // [Archive_X, Archive_F, Archive_member_no]=UpdateArchive(Archive_X, Archive_F, X, Particles_F, Archive_member_no);
            Object tmp[] = new Object[3];
            tmp = UpdateArchive(Archive_X, Archive_F, X, Particles_F, Archive_member_no);
            Archive_X = (double[][]) tmp[0];
            Archive_F = (double[][]) tmp[1];
            Archive_member_no = (int) tmp[2];

            double Archive_mem_ranks[] = null;

            // if Archive_member_no>ArchiveMaxSize
            if (Archive_member_no > ArchiveMaxSize) {
                // Archive_mem_ranks=RankingProcess(Archive_F, ArchiveMaxSize, obj_no);
/* bu blok dogrulanmadi */
                Archive_mem_ranks = RankingProcess(Archive_F, ArchiveMaxSize, obj_no);
                Object tmp2[] = HandleFullArchive(Archive_X, Archive_F, Archive_member_no, Archive_mem_ranks, ArchiveMaxSize);
                Archive_X = (double[][]) tmp2[0];
                Archive_F = (double[][]) tmp2[1];
                Archive_mem_ranks = (double[]) tmp2[2];
                Archive_member_no = (int) tmp2[3];

            } else {
                // Archive_mem_ranks=RankingProcess(Archive_F, ArchiveMaxSize, obj_no);
                Archive_mem_ranks = RankingProcess(Archive_F, ArchiveMaxSize, obj_no);
            }
            // Archive_mem_ranks
            Archive_mem_ranks = RankingProcess(Archive_F, ArchiveMaxSize, obj_no);
            // index=RouletteWheelSelection(1./Archive_mem_ranks);
            int index = RouletteWheelSelection(elementWiseDivision(onesRowVector(getSizeOfDim(Archive_mem_ranks, 2)), Archive_mem_ranks));
            // if index==-1
            if (index == -1) {
                // index=1;
                index = 1;
            }
            // Food_fitness=Archive_F(index,:);
            Food_fitness = getIthRow(Archive_F, index);
            // Food_pos=Archive_X(index,:)';
            Food_pos = transposeRow(getIthRow(Archive_X, index));
            // index=RouletteWheelSelection(Archive_mem_ranks);
            index = RouletteWheelSelection(Archive_mem_ranks);
            if (index == -1) {
                // index=1;
                index = 1;
            }
            // Enemy_fitness=Archive_F(index,:);
            Enemy_fitness = getIthRow(Archive_F, index);
            // Enemy_pos=Archive_X(index,:)';
            Enemy_pos = transposeRow(getIthRow(Archive_X, index));
            // neighbours_no
            int neighbours_no;
            double Neighbours_V[][] = null;
            double Neighbours_X[][] = null;
            // for i=1:N
            for (int i = 0; i < N; i++) {
                // index=0;
                index = 0;
                // neighbours_no=0;
                neighbours_no = 0;
                // clear Neighbours_V
                Neighbours_V = null;
                // clear Neighbours_X
                Neighbours_X = null;
                double Dist[] = null;
                // for j=1:N
                for (int j = 0; j < N; j++) {
                    // Dist=distance(X(:,i),X(:,j));
                    Dist = distance(getJthColumn(X, i), getJthColumn(X, j));
                    // if (all(Dist<=r) && all(Dist~=0))
                    if (isAllNonZero(isLessThenAndEqual(Dist, r[0]))
                            && isAllNonZero(isNotEqual(Dist, zerosRowVector(Dist.length)))) {
                        // index=index+1;
                        index++;
                        // neighbours_no=neighbours_no+1;
                        neighbours_no++;
                        // Neighbours_V(:,index)=DeltaX(:,j);
                        setColumn(Neighbours_V, index, getJthColumn(DeltaX, j));
                        // Neighbours_X(:,index)=X(:,j);
                        setColumn(Neighbours_X, index, getJthColumn(X, j));
                    }
                }
                // S=zeros(dim,1);
                double S[][] = zerosColumnVector(dim);
                // if neighbours_no>1
                if (neighbours_no > 1) {
                    // for k=1:neighbours_no
                    for (int j = 0; j < neighbours_no; j++) {
                        // S=S+(Neighbours_X(:,k)-X(:,i));
                        S = elementWiseAddition(S,
                                elementWiseSubstraction(
                                        getJthColumn(Neighbours_X, j),
                                        getJthColumn(X, i)));
                    }
                    // S=-S;
                    S = elementWiseProduct(S, -1);
                } else {
                    // S=zeros(dim,1);
                    S = zerosColumnVector(dim);
                }
                double A[][];
                // if neighbours_no>1
                if (neighbours_no > 1) {
                    // A=(sum(Neighbours_V')')/neighbours_no;
                    A = elementWiseDivision(
                            transposeRow(sum(transposeMatrix(Neighbours_V))), neighbours_no);
                } else {
                    // A=DeltaX(:,i);
                    A = getJthColumn(DeltaX, i);
                }
                double C_temp[][] = null;
                // if neighbours_no>1
                if (neighbours_no > 1) {
                    // C_temp=(sum(Neighbours_X')')/neighbours_no;
                    C_temp = elementWiseDivision(
                            transposeRow(sum(transposeMatrix(Neighbours_X))), neighbours_no);
                } else {
                    // C_temp=X(:,i);
                    C_temp = getJthColumn(X, i);
                }
                // C=C_temp-X(:,i);
                double C[][] = elementWiseSubstraction(C_temp, getJthColumn(X, i));
                // Dist2Attraction=distance(X(:,i),Food_pos(:,1));
                double Dist2Attraction[] = distance(getJthColumn(X, i), getJthColumn(Food_pos, 1));
                double F[][] = null;
                // if all(Dist2Attraction<=r)
                if (isAllNonZero(isLessThenAndEqual(Dist2Attraction, r[0]))) {
                    // F=Food_pos-X(:,i);
                    F = elementWiseSubstraction(Food_pos, getJthColumn(X, i));
                    // iter;
                    System.out.print(iter);
                } else {
                    // F=0;
                    F = zerosColumnVector(dim);
                }
                double E[][] = null;
                // Dist=distance(X(:,i),Enemy_pos(:,1));
                Dist = distance(getJthColumn(X, i), getJthColumn(Enemy_pos, 1));
                // if all(Dist<=r)
                if (isAllNonZero(isLessThenAndEqual(Dist, r[0]))) {
                    // E=Enemy_pos+X(:,i);
                    E = getJthColumn(Enemy_pos, i);
                } else {
                    // E=zeros(dim,1);
                    E = zerosMatrix(dim, 1);
                }
                // for tt=1:dim
                for (int tt = 0; tt < dim; tt++) {
                    // if X(tt,i)>ub(tt)
                    if (X[tt][i] > ub[tt][0]) {
                        // X(tt,i)=lb(tt);
                        X[tt][i] = lb[tt][0];
                        // DeltaX(tt,i)=rand;
                        DeltaX[tt][i] = Math.random();
                    }
                    // if X(tt,i)<lb(tt)
                    if (X[tt][i] < lb[tt][0]) {
                        // X(tt,i)=ub(tt);
                        X[tt][i] = ub[tt][0];
                        // DeltaX(tt,i)=rand;
                        DeltaX[tt][i] = Math.random();
                    }
                }
                if (isAnyNonZero(isGreaterThen(Dist2Attraction, r[0]))) {
                    if (neighbours_no > 1) {
                        for (int j = 0; j < dim; j++) {
                            DeltaX[j][i] = w * DeltaX[j][i]
                                    + Math.random() * A[j][1]
                                    + Math.random() * C[j][1]
                                    + Math.random() * S[j][1];
                            if (DeltaX[j][i] > V_max) {
                                DeltaX[j][i] = V_max;
                            }
                            if (DeltaX[j][i] < -V_max) {
                                DeltaX[j][i] = -V_max;
                            }
                            X[j][i] = X[j][i] + DeltaX[j][i];
                        }
                    } else {
                        setColumn(X, i,
                                elementWiseAddition(getJthColumn(X, i),
                                        elementWiseProduct(transposeRow(Levy(dim)), getJthColumn(X, i))));
                        setColumn(DeltaX, i, zerosColumnVector(DeltaX[0].length));

                    }
                } else {
                    for (int j = 0; j < dim; j++) {
                        DeltaX[j][i] = s * S[j][1]
                                + a * A[j][1]
                                + c * C[j][1]
                                + f * F[j][1]
                                + e * E[j][1]
                                + w * DeltaX[j][i];
                        if (DeltaX[j][i] > V_max) {
                            DeltaX[j][i] = V_max;
                        }
                        if (DeltaX[j][i] < -V_max) {
                            DeltaX[j][i] = -V_max;
                        }
                        X[j][1] = X[j][i] + DeltaX[j][i];
                    }
                }
                double Flag4ub[][] = isGreaterThen(getJthColumn(X, i), transposeMatrix(ub));
                double Flag4lb[][] = isLessThen(getJthColumn(X, i), transposeMatrix(lb));
                setColumn(X, i,
                        elementWiseAddition(elementWiseAddition(
                                elementWiseProduct(getJthColumn(X, i), logicNOT(logicOR(Flag4ub, Flag4lb))),
                                elementWiseProduct(transposeMatrix(ub), Flag4ub)),
                                elementWiseProduct(transposeMatrix(lb), Flag4lb)));
                

            }
        }
        System.out.println("At the iteration " + iter + " there are " + Archive_member_no + "  non-dominated solutions in the archive" );
    }

    public static double[][] logicNOT(double m1[][]) {
        double res[][] = new double[m1.length][m1[0].length];
        for (int i = 0; i < m1.length; i++) {
            for (int j = 0; j < m1[0].length; j++) {
                res[i][j] = (m1[i][j] == 0) ? 1 : 0;
            }
        }
        return res;
    }

    public static double[][] logicOR(double m1[][], double m2[][]) {
        double res[][] = new double[m1.length][m1[0].length];
        for (int i = 0; i < m1.length; i++) {
            for (int j = 0; j < m1[0].length; j++) {
                res[i][j] = (m1[i][j] == 0 && m2[i][j] == 0) ? 0 : 1;
            }
        }
        return res;
    }

    public static double[] distance(double m1[][], double m2[][]) {
        double res[] = new double[m1.length];
        for (int i = 0; i < m1.length; i++) {
            // o(1,i)=sqrt((a(i)-b(i))^2);
            if (getSizeOfDim(m1, 2) == 1) {
                res[i] = m1[i][0] - m2[i][0];
            } else if (getSizeOfDim(m1, 2) != 1) {
                /* burasi dolacak */
            }

        }
        return res;
    }

    public static Object[] HandleFullArchive(double Archive_X[][], double Archive_F[][], int Archive_member_no, double Archive_mem_ranks[], int ArchiveMaxSize) {
        Object results[] = new Object[4];
        double Archive_X_Chopped[][] = Archive_X;
        double Archive_F_Chopped[][] = Archive_F;
        double Archive_mem_ranks_updated[] = Archive_mem_ranks;

        for (int i = 0; i < getSizeOfDim(Archive_F, 1) - ArchiveMaxSize; i++) {
            int index = RouletteWheelSelection(Archive_mem_ranks);

            Archive_X_Chopped = mergeRows(getRows(Archive_X_Chopped, 0, index - 1), getRows(Archive_X_Chopped, index + 1, Archive_member_no));
            Archive_F_Chopped = mergeRows(getRows(Archive_F_Chopped, 0, index - 1), getRows(Archive_F_Chopped, index + 1, Archive_member_no));
            Archive_mem_ranks_updated = mergeColumns(getSubRow(Archive_mem_ranks, 0, index - 1), getSubRow(Archive_mem_ranks, index + 1, Archive_member_no));
            Archive_member_no--;
        }
        results[0] = Archive_X_Chopped;
        results[1] = Archive_F_Chopped;
        results[2] = Archive_mem_ranks_updated;
        results[3] = Archive_member_no;

        return results;
    }

    public static double[] RankingProcess(double Archive_F[][], int ArchiveMaxSize, int obj_no) {

        double my_min[] = null;
        double my_max[] = null;

        // my_min=min(Archive_F);
        my_min = min(Archive_F);
        // my_max=max(Archive_F);
        my_max = max(Archive_F);
        // r=(my_max-my_min)/(20);
        double r[] = elementWiseDivision(elementWiseSubstraction(my_max, my_min), new double[]{20, 20});
        // ranks=zeros(1,size(Archive_F,1));
        double ranks[] = zerosRowVector(getSizeOfDim(Archive_F, 1));
        // for i=1:size(Archive_F,1)
        for (int i = 0; i < ranks.length; i++) {
            // ranks(i)=0;
            ranks[0] = 0;
            // for j=1:size(Archive_F,1)
            for (int j = 0; j < ranks.length; j++) {
                // flag=0;
                int flag = 0;
                // for k=1:obj_no
                for (int k = 0; k < obj_no; k++) {
                    // if (abs(Archive_F(j,k)-Archive_F(i,k))<r(k))
                    if (Math.abs(Archive_F[j][k] - Archive_F[i][k]) < r[k]) {
                        // flag=flag+1;
                        flag++;
                    }
                }
                // if flag==obj_no
                if (flag == obj_no) {
                    // ranks(i)=ranks(i)+1;
                    ranks[i] = ranks[i] + 1;
                }
            }
        }
        return ranks;
    }

    public static double[] max(double m[][]) {
        double res[] = new double[2];
        res = m[0];
        for (int i = 0; i < m.length; i++) {
            if (m[i][0] > res[0]) {
                res[0] = m[i][0];
            }
            if (m[i][1] > res[1]) {
                res[1] = m[i][1];
            }
        }
        return res;
    }

    public static double[] min(double m[][]) {
        double res[] = new double[2];
        res = m[0];
        for (int i = 0; i < m.length; i++) {
            if (m[i][0] < res[0]) {
                res[0] = m[i][0];
            }
            if (m[i][1] < res[1]) {
                res[1] = m[i][1];
            }
        }
        return res;
    }

    public static int RouletteWheelSelection(double[] weights) {
        double[] accumulation = new double[weights.length];
        int chosen_index = 0;
        double kumulatifToplam = 0;
        double p;

        for (int i = 0; i < weights.length; i++) {
            kumulatifToplam = kumulatifToplam + weights[i];
            accumulation[i] = kumulatifToplam;
        }
        p = Math.random() * accumulation[accumulation.length];
        for (int index = 0; index < accumulation.length; index++) {
            if (accumulation[index] > p) {
                chosen_index = index;
                break;
            }
        }

        return chosen_index;
    }

    // dogrulandi
    public static double[] Levy(int dim) {
        // beta=3/2;
        double beta = 1.5;
        // sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2))) ^(1/beta);
        double sigma = Math.pow(gamma(1 + beta) * Math.sin(3.1415 * beta / 2) / (gamma((1 + beta) / 2) * beta * Math.pow(2, (beta - 1) / 2)), 1 / beta);
        // u=randn(1,d)*sigma;
        double u[] = elementWiseProduct(randn(1, dim), sigma)[0];
        // v=randn(1,d);
        double v[] = randn(1, dim)[0];
        // step=u./abs(v).^(1/beta);
        double step[] = elementWiseDivision(u, elementWisePower(elementWiseAbs(v), 1 / beta));
        // o=0.01*step;
        return elementWiseProduct(step, 0.01);
    }

    // dogrulandi
    public static double logGamma(double x) {
        double tmp = (x - 0.5) * Math.log(x + 4.5) - (x + 4.5);
        double ser = 1.0 + 76.18009173 / (x + 0) - 86.50532033 / (x + 1)
                + 24.01409822 / (x + 2) - 1.231739516 / (x + 3)
                + 0.00120858003 / (x + 4) - 0.00000536382 / (x + 5);
        return tmp + Math.log(ser * Math.sqrt(2 * Math.PI));
    }

    // dogrulandi
    public static double gamma(double x) {
        return Math.exp(logGamma(x));
    }

    // dogrulandi
    public static double[][] randn(int x, int y) {
        double rands[][] = new double[x][y];
        for (int i = 0; i < x; i++) {
            for (int j = 0; j < y; j++) {
                rands[i][j] = Math.pow(-1, (int) (20 * Math.random())) * Math.random();
            }
        }
        return rands;
    }

    // dogrulandi
    public static Object[] UpdateArchive(double Archive_X[][], double Archive_F[][], double Particles_X[][], double Particles_F[][], int Archive_member_no) {
        Object results[] = new Object[3];
        // Archive_X_temp=[Archive_X ; Particles_X'];
        double Archive_X_temp[][] = mergeRows(Archive_X, transposeMatrix(Particles_X));
        // Archive_F_temp=[Archive_F ; Particles_F];
        double Archive_F_temp[][] = mergeRows(Archive_F, Particles_F);
        // o=zeros(1,size(Archive_F_temp,1));
        double o[] = zerosRowVector(getSizeOfDim(Archive_F_temp, 1));
        // for i=1:size(Archive_F_temp,1)
        for (int i = 0; i < getSizeOfDim(Archive_F_temp, 1); i++) {
            // o(i)=0;
            o[i] = 0;
            // for j=1:i-1
            for (int j = 0; j < i; j++) {
                // if any(Archive_F_temp(i,:) ~= Archive_F_temp(j,:))
                if (isAnyNonZero(isNotEqual(getIthRow(Archive_F_temp, i), getIthRow(Archive_F_temp, j)))) {
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
        double Archive_X_updated[][] = null;
        double Archive_F_updated[][] = null;

        // for i=1:size(Archive_X_temp,1)
        for (int i = 0; i < getSizeOfDim(Archive_X_temp, 1); i++) {
            // if o(i)==0
            if (o[i] == 0) {
                // Archive_X_updated(Archive_member_no,:)=Archive_X_temp(i,:);
                Archive_X_updated = mergeRow(Archive_X_updated, getIthRow(Archive_X_temp, i));
                // Archive_F_updated(Archive_member_no,:)=Archive_F_temp(i,:);
                Archive_F_updated = mergeRow(Archive_F_updated, getIthRow(Archive_F_temp, i));
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

    // dogrulandi
    public static boolean dominates(double x[], double y[]) {
        return isAllNonZero(isLessThenAndEqual(x, y)) && isAnyNonZero(isLessThen(x, y));

    }

    // dogrulandi
    public static double[] ZDT1(double x[]) {
        // o = [0, 0];
        double o[] = zerosRowVector(2);
        // dim = length(x);
        int dim = x.length;
        // g = 1 + 9*sum(x(2:dim))/(dim-1);
        double g = 1 + (9 * sum(getSubRow(x, 1, dim - 1)) / (dim - 1));
        o[0] = x[0];
        o[1] = g * (1 - (Math.sqrt(x[0] / g)));
        return o;
    }

    // dogrulandi
    public static double sum(double v[]) {
        double res = 0.0;
        for (int i = 0; i < v.length; i++) {
            res += v[i];
        }
        return res;
    }

    // dogrulandi
    public static double[] sum(double m[][]) {
        double res[] = new double[m[0].length];
        for (int j = 0; j < m[0].length; j++) {
            for (int i = 0; i < m.length; i++) {
                res[j] += m[j][i];
            }
        }
        return res;
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
        return transposeMatrix(X);
    }

    // dogrulandi
    public static double[][] transposeMatrix(double[][] m) {
        double[][] temp = new double[m[0].length][m.length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                temp[j][i] = m[i][j];
            }
        }
        return temp;
    }

    // dogrulandi
    public static double[] transposeColumn(double m[][]) {
        double arr[] = new double[m.length];
        for (int i = 0; i < m.length; i++) {
            arr[i] = m[i][0];
        }
        return arr;
    }

    // dogrulandi
    public static double[][] transposeRow(double vec[]) {
        double m[][] = new double[vec.length][1];
        for (int i = 0; i < vec.length; i++) {
            m[i][0] = vec[i];
        }
        return m;
    }

    public static void setRow(double m[][], int row, double v[]) {
        for (int i = 0; i < m[row].length; i++) {
            m[row][i] = v[i];
        }
    }

    public static void setColumn(double m[][], int column, double v[][]) {
        for (int i = 0; i < m.length; i++) {
            m[i][column] = v[i][0];
        }
    }

    // dogrulandi
    public static double[] getIthRow(double matrix[][], int I) {
        double arr[] = new double[matrix[0].length];
        for (int i = 0; i < matrix[0].length; i++) {
            arr[i] = matrix[I][i];
        }
        return arr;
    }

    // dogrulandi
    public static double[] getSubRow(double row[], int start, int end) {
        double sub[] = new double[end - start + 1];
        for (int i = 0, j = start; i < sub.length; i++, j++) {
            sub[i] = row[j];
        }
        return sub;
    }

    public static double[][] getSubCloumn(double column[][], int start, int end) {
        double sub[][] = new double[end - start + 1][1];
        for (int i = start; i < sub.length; i++) {
            sub[i][0] = column[i][0];
        }
        return sub;
    }

    public static double[][] getRows(double matrix[][], int rowStart, int rowEnd) {
        double res[][] = new double[rowEnd - rowStart + 1][matrix[0].length];
        for (int i = 0, j = rowStart; i < res.length; i++, j++) {
            for (int k = 0; k < matrix[0].length; k++) {
                res[i][k] = matrix[j][k];
            }
        }
        return res;
    }

    // dogrulandi
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

    // dogrulandi
    public static double[][] mergeRows(double m1[][], double m2[][]) {
        double m3[][] = new double[m1.length + m2.length][m1[0].length];
        for (int i = 0; i < m1.length; i++) {
            for (int j = 0; j < m1[0].length; j++) {
                m3[i][j] = m1[i][j];
            }
        }
        for (int i = 0; i < m2.length; i++) {
            for (int j = 0; j < m2[0].length; j++) {
                m3[i + m1.length][j] = m2[i][j];
            }
        }
        return m3;
    }

    public static double[] mergeColumns(double m1[], double m2[]) {
        double res[] = new double[m1.length + m2.length];
        for (int i = 0; i < m1.length; i++) {
            res[i] = m1[i];
        }
        for (int i = 0; i < m2.length; i++) {
            res[i + m1.length] = m2[i];
        }
        return res;
    }

    public static double[][] mergeColumns(double m1[][], double m2[][]) {
        double res[][] = null;
        if (m1 == null) {
            res = new double[m2.length][m2[0].length];
            res = m2;
            return res;
        } else {
            res = new double[m2.length][m2[0].length + m1[0].length];
            for (int i = 0; i < m1.length; i++) {
                for (int j = 0; j < m1[0].length; j++) {
                    res[i][j] = m1[i][j];
                }
            }
            for (int i = 0; i < m2.length; i++) {
                for (int j = 0; j < m2[0].length; j++) {
                    res[i][m1[0].length + j] = m2[i][j];
                }
            }
            return res;
        }
    }

    // dogrulandi
    public static double[][] mergeRow(double m1[][], double m2[]) {
        int cn = m2.length;
        int rn = 0;
        if (m1 != null) {
            rn = m1.length;
            cn = m1[0].length;
        }
        double m3[][] = new double[rn + 1][cn];
        if (m1 != null) {
            for (int i = 0; i < rn; i++) {
                m3[i] = m1[i];
            }
        }
        m3[rn] = m2;
        return m3;
    }

    public static double[][] elementWiseDivision(double m[][], double scalar) {
        double res[][] = new double[m.length][m[0].length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                res[i][j] = m[i][j] / scalar;
            }
        }
        return res;
    }

    public static double[] elementWiseDivision(double m[], double m2[]) {
        double res[] = new double[m.length];
        for (int i = 0; i < m.length; i++) {
            res[i] = m[i] / m2[i];
        }
        return res;
    }

    public static double[][] elementWiseProduct(double m[][], double scalar) {
        double res[][] = new double[m.length][m[0].length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                res[i][j] = m[i][j] * scalar;
            }
        }
        return res;
    }

    public static double[][] elementWiseProduct(double m1[][], double m2[][]) {
        double res[][] = new double[m1.length][m1[0].length];
        for (int i = 0; i < m1.length; i++) {
            for (int j = 0; j < m1[0].length; j++) {
                res[i][j] = m1[i][j] * m2[i][j];
            }
        }
        return res;
    }

    public static double[] elementWiseProduct(double vec[], double scalar) {
        double res[] = new double[vec.length];
        for (int i = 0; i < vec.length; i++) {
            res[i] = vec[i] * scalar;
        }
        return res;
    }

    public static double[][] elementWiseAddition(double m[][], double scalar) {
        double res[][] = new double[m.length][m[0].length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                res[i][j] = m[i][j] + scalar;
            }
        }
        return res;
    }

    public static double[][] elementWiseAddition(double m[][], double n[][]) {
        double res[][] = new double[m.length][m[0].length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                res[i][j] = m[i][j] + n[i][j];
            }
        }
        return res;
    }

    public static double[][] elementWiseSubstraction(double m[][], double scalar) {
        double res[][] = new double[m.length][m[0].length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                res[i][j] = m[i][j] - scalar;
            }
        }
        return res;
    }

    public static double[][] elementWiseSubstraction(double m[][], double n[][]) {
        double res[][] = new double[m.length][m[0].length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                res[i][j] = m[i][j] - n[i][j];
            }
        }
        return res;
    }

    public static double[] elementWiseSubstraction(double m[], double n[]) {
        double res[] = new double[m.length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m.length; j++) {
                res[i] = m[i] - n[i];
            }
        }
        return res;
    }

    public static double[] elementWisePower(double m[], double m2) {
        double res[] = new double[m.length];
        for (int i = 0; i < m.length; i++) {
            res[i] = Math.pow(m[i], m2);
        }
        return res;
    }

    public static double[] elementWisePower(double m1[], double m2[]) {
        double res[] = new double[m1.length];
        for (int i = 0; i < res.length; i++) {
            res[i] = Math.pow(m1[i], m2[i]);
        }
        return res;
    }

    public static double[] elementWiseAbs(double m1[]) {
        double res[] = new double[m1.length];
        for (int i = 0; i < res.length; i++) {
            res[i] = Math.abs(m1[i]);
        }
        return res;
    }

    // dogrulandi
    public static int getSizeOfDim(double m[], int dim) {
        int size = 0;
        switch (dim) {
            case 1:
                size = m.length;
                break;
        }
        return size;
    }

    // dogrulandi
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

    // dogrulandi
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

    public static boolean[] isAllNonZero(double m[][]) {
        boolean result[] = new boolean[m.length];
        for (int i = 0; i < m[0].length; i++) {
            result[i] = isAllNonZero(m, i);
        }
        return result;
    }

    // dogrulandi
    public static boolean isAllNonZero(double m[]) {
        for (int i = 0; i < m.length; i++) {
            if (m[i] == 0) {
                return false;
            }
        }
        return true;
    }

    // dogrulandi
    public static boolean isAllNonZero(double m[][], int c) {
        for (int i = 0; i < m.length; i++) {
            if (m[i][c] == 0) {
                return false;
            }
        }
        return true;
    }

    public static boolean[] isAnyNonZero(double m[][]) {
        boolean result[] = new boolean[m.length];
        for (int i = 0; i < m[0].length; i++) {
            result[i] = isAnyNonZero(m, i);
        }
        return result;
    }

    public static boolean isAnyNonZero(double m[]) {
        for (int i = 0; i < m.length; i++) {
            if (m[i] != 0) {
                return true;
            }
        }
        return false;
    }

    public static boolean isAnyNonZero(double m[][], int c) {
        for (int i = 0; i < m.length; i++) {
            if (m[i][c] != 0) {
                return true;
            }
        }
        return false;
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

    // dogrulandi
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

    // dogrulandi
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

    public static double[][] isLessThen(double m1[][], double m2[][]) {
        double res[][] = new double[m1.length][1];
        for (int i = 0; i < m1.length; i++) {
            if (m1[i][0] < m2[i][0]) {
                res[i][0] = 1;
            } else {
                res[i][0] = 0;
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

    // dogrulandi
    public static double[] isGreaterThen(double m1[], double m2[]) {
        double res[] = new double[m1.length];
        for (int i = 0; i < m1.length; i++) {
            if (m1[i] > m2[i]) {
                res[i] = 1;
            } else {
                res[i] = 0;
            }
        }
        return res;
    }

    public static double[][] isGreaterThen(double m1[][], double m2[][]) {
        double res[][] = new double[m1.length][1];
        for (int i = 0; i < m1.length; i++) {
            if (m1[i][0] > m2[i][0]) {
                res[i][0] = 1;
            } else {
                res[i][0] = 0;
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

    // dogrulandi
    public static void print(double m[]) {
        System.out.print("[\t");
        for (int i = 0; i < m.length; i++) {
            System.out.print(m[i]);
            System.out.print("\t");
        }
        System.out.print("\t]");
    }

    public static void print(boolean m[]) {
        System.out.print("[\t");
        for (int i = 0; i < m.length; i++) {
            System.out.print(m[i]);
            System.out.print("\t");
        }
        System.out.print("\t]");
    }

    // dogrulandi
    public static void print(double m[][]) {

        for (int i = 0; i < m.length; i++) {
            System.out.print("[\t");
            for (int j = 0; j < m[0].length; j++) {
                System.out.print(m[i][j]);
                System.out.print("\t");
            }
            System.out.print("\t]\n");
        }

    }

}
