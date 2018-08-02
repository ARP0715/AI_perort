//
// Created by Admin on 2018/07/25.
//

/********************************************************************
遺伝的アルゴリズム
********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <math.h>

#define POP_SIZE 5   /* 個体数 (必ず奇数に設定) */
#define G_LENGTH 10  /* 個体の遺伝子型のビット数 */
#define MAX_GEN  20  /* 世代数 */
#define M_RATE   0.1 /* 突然変異率 (0～1) */

double avg_fit;  /* 平均適応度  */
double max_fit;
int t;
FILE *fp;


/********************************************************************
乱数の発生 (Seedの決定)
********************************************************************/
void init_rnd()
{
    srand((unsigned int)time(NULL));
}

/********************************************************************
乱数の発生 (0～1の乱数)
********************************************************************/
double Random()
{
    return((double)rand() / RAND_MAX);
}

/********************************************************************
遺伝子の初期化
引数 gene[p][i] : 遺伝子pのi番目の成分
********************************************************************/
void init_gene(int gene[POP_SIZE][G_LENGTH])
{
    int p, i;

    /* 乱数の初期化 */
    init_rnd();

    /* 遺伝子を初期化  0～1の乱数を発生し、0.5以上なら1
    0.5未満なら0 */
    printf("<< 初期個体群 >>\n");

    /*** この部分を自分で書く  ***/
    for (p = 0; p < POP_SIZE; p++) {//ループ処理前に0をいれ、POP_SIZEの処理を繰り返す。
        for (i = 0; i < G_LENGTH; i++) {//ループ処理前に0をいれ、G_LENGTHの処理を繰り返す。

            if (Random() < 0.5) gene[p][i] = 0; //もしRandom()が0.5より小さければ0となる。
            else gene[p][i] = 1; //それ以外は1となる。
        }
    }
}

/********************************************************************
適応度の計算
引数 gene[p][i] : 遺伝子pのi番目の成分
fitness[p] : 遺伝子pの適応度
********************************************************************/
void calc_fitness(int gene[POP_SIZE][G_LENGTH], double fitness[POP_SIZE])
{
    int p, i;

    /* 適応度の計算 前半の5bitは0の数 後半の5bitは1の数 */

    /*** この部分を自分で書く  ***/
    for (p = 0; p < POP_SIZE; p++) {//ループ処理前に0をいれ、POP_SIZEの処理を繰り返す。

        fitness[p] = 0;//fitnessに0をいれる。
        for (i = 0; i < 5; i++) {//iは0からはじまり５になるまで繰り返す。
            if (gene[p][i] == 0) fitness[p] += 1;// gene[p][i] == 0ならfitness[p]に1を足す。
        }
        for (; i < G_LENGTH; i++) {//iをG_LENGTHの処理を繰り返す。
            if (gene[p][i] == 1) fitness[p] += 1;//もしgene[p][i] == 1ならfitness[p]に1を足す。
        }
    }
}

/**********************************************************************
遺伝子の表示 & 最大適応度・平均適応度の計算 & ファイルへの書き出し
引数 t          : 世代数
gene[p][i] : 遺伝子pのi番目の成分
fitness[p] : 遺伝子pの適応度
*fp        : ファイルポインタ
**********************************************************************/
void show_gene(int t, int gene[POP_SIZE][G_LENGTH], double fitness[POP_SIZE], FILE *fp) {
    //int t, int gene[POP_SIZE][G_LENGTH], double fitness[POP_SIZE], FILE *fp
    int p, i;
    double avg_fit;  /* 平均適応度  */
    double max_fit;  /* 最大適応度  */

    /* 個体の値、適応度の表示 */

    /*** この部分を自分で書く  ***/
    avg_fit = 0.0;//avg_fitに0.0をいれる。
    max_fit = 0.0;//max_fitに0.0をいれる。
    /* 平均・最大適応度の計算 */
    for (p = 0; p < POP_SIZE; p++) {//ループ処理前に0をいれ、POP_SIZEの処理を繰り返す。
        if (max_fit < fitness[p]) max_fit = fitness[p];//もしmax_fitはfitness[p]より小さければ、max_fitにfitness[p]をコピーする。
        avg_fit += fitness[p]; // fitness[p]からavg_fitを足す。
    }
    avg_fit /= POP_SIZE;//POP_SIZEからavg_fitを割る。



    /* 平均・最大適応度の表示 */
    printf("平均適応度 : %lf\n", avg_fit);
    printf("最大適応度 : %lf\n", max_fit);

    /* 平均・最大適応度をファイルに書き込む */
    fprintf(fp, "%d %lf %lf\n", t, avg_fit, max_fit);

}

/**********************************************************************
個体番号 p1 と p2 の適応度と遺伝子を交換
引数 p1, p2     : 遺伝子の番号
gene[p][i] : 遺伝子pのi番目の成分
fitness[p] : 遺伝子pの適応度
**********************************************************************/
void swap_gene(int p1, int p2, int gene[POP_SIZE][G_LENGTH],
               double fitness[POP_SIZE])
{
    int tmp;
    double f;
    int i;

    /* 遺伝子型の交換 (遺伝子p1と遺伝子p2の値を入れ替える) */

    /*** この部分を自分で書く  ***/
    for (i = 0; i < G_LENGTH; i++) { //ループ処理前に0をいれ、G_LENGTHの処理を繰り返す。
        tmp = gene[p1][i];//gene[p1][i]をtmpにコピーする。
        gene[p1][i] = gene[p2][i];//gene[p2][i]をgene[p1][i]にコピーする。
        gene[p2][i] = tmp; //tmpをgene[p2][i]にコピーする。
    }

    /* 適応度の交換 (遺伝子p1と遺伝子p2の適応度の値を入れ替える) */

    /*** この部分を自分で書く  ***/
    f = fitness[p1];//fにfitness[p1]をいれる。
    fitness[p1] = fitness[p2];//fitness[p2]をfitness[p1]にコピーする。
    fitness[p2] = f;//fをfitness[p2]にコピーする。


}

/**********************************************************************
個体番号 p1 の適応度と遺伝子型を p2 にコピー
引数 p1, p2     : 遺伝子の番号
gene[p][i] : 遺伝子pのi番目の成分
fitness[p] : 遺伝子pの適応度
**********************************************************************/
void copy_gene(int p1, int p2, int gene[POP_SIZE][G_LENGTH],
               double fitness[POP_SIZE])
{
    int i;

    /* 遺伝子のコピー (遺伝子p1を遺伝子p2にコピーする) */

    /*** この部分を自分で書く  ***/
    for (i = 0; i < G_LENGTH; i++) {//ループ処理前に0をいれ、G_LENGTHの処理を繰り返す。
        gene[p2][i] = gene[p1][i]; //gene[p1][i]をgene[p2][i]にコピーする。
    }

    /* 適応度のコピー (遺伝子p1の適応度を遺伝子p2の適応度にコピーする)*/

    /*** この部分を自分で書く  ***/
    fitness[p2] = fitness[p1];//fitness[p1]をfitness[p2]にコピーする。

}

/**********************************************************************
エリート保存
(最小適応度の個体に最大適応度の個体のデータをコピー)
引数 gene[p][i] : 遺伝子pのi番目の成分
fitness[p] : 遺伝子pの適応度
**********************************************************************/
void elite(int gene[POP_SIZE][G_LENGTH], double fitness[POP_SIZE])
{
    int p, i;
    double max_fitness = fitness[0];
    double min_fitness = fitness[0];
    int max_p = 0;
    int min_p = 0;

    /* 最大適応度の個体(max_p)と最小適応度の個体(min_p)を見つける */

    /*** この部分を自分で書く  ***/
    for (p = 1; p < POP_SIZE; p++) {//ループ処理前に1をいれ、POP_SIZEの処理を繰り返す。
        if (max_fitness < fitness[p]) {  //もしmax_fitnessがfitness[p]よりちいさければ、
            max_fitness = fitness[p]; //fitness[p]をmax_fitnessにコピーする。
            max_p = p;//pをmax_pにコピーする。
        }
        if (min_fitness > fitness[p]) {//もしmin_fitnessがfitness[p]より大きければ、
            min_fitness = fitness[p];//fitness[p]をmin_fitnessにコピーする。
            min_p = p;//pをmin_pにコピーする。
        }
    }



    /* 最小適応度の個体に最大適応度の個体をコピー */
    copy_gene(max_p, min_p, gene, fitness);
    /* 最大適応度の個体を0番目に移動 */
    swap_gene(0, max_p, gene, fitness);
}

/**********************************************************************
ルーレット選択
引数 gene[p][i] : 遺伝子pのi番目の成分
fitness[p] : 遺伝子pの適応度
**********************************************************************/
void reproduction(int gene[POP_SIZE][G_LENGTH], double fitness[POP_SIZE])
{
    double sum_of_fitness; /* 個体の適応度の総和 */
    double border;         /* ルーレット上の個体間の境界 */
    double r;              /* ルーレット上の選択位置 */
    int p, i;               /* 選ばれた個体の番号 */
    int num;               /* 0 <= num <= POP_SIZE-1 */
    int new_gene[POP_SIZE][G_LENGTH];

    /* ルーレットの1周分 sum_of_fitness を求める */
    sum_of_fitness = 0.0;
    for (p = 0; p<POP_SIZE; p++) {
        sum_of_fitness += fitness[p];
    }

    /* ルーレットを POP_SIZE 回だけ回して次世代の個体を選ぶ */
    for (p = 1; p<POP_SIZE; p++) {
        /* ルーレットを回して場所を選ぶ
        r : 選ばれた位置 (0 <= r <= sum_of_fitness) */
        r = sum_of_fitness*Random();
        /* 選ばれた場所に該当する個体が何番か調べる
        num : 選ばれた個体の番号 (0 <= num <= POP_SIZE-1) */
        num = 0;
        border = fitness[0]; /* 個体間の境界 */
        while (border<r) {
            num++;
            border += fitness[num];
        }

        /* 遺伝子の代入 */
        for (i = 0; i<G_LENGTH; i++) {
            new_gene[p][i] = gene[num][i];
        }
    }

    /* 遺伝子のコピー */
    for (p = 1; p<POP_SIZE; p++) {
        for (i = 0; i<G_LENGTH; i++) {
            gene[p][i] = new_gene[p][i];
        }
    }
}

/**********************************************************************
一点交叉
引数 gene[p][i] : 遺伝子pのi番目の成分
**********************************************************************/
void crossover(int gene[POP_SIZE][G_LENGTH])
{
    int gene1[G_LENGTH]; /* 親1の遺伝子型 */
    int gene2[G_LENGTH]; /* 親2の遺伝子型 */
    int i, j;
    int c_pos;   /* 交叉位置 (1 <= c_pos <= G_LENGTH-1) */

    /* 交叉位置を1～G_LENGTH-1の範囲でランダムに決め、
    それより後ろを入れ替える。
    gene[1]とgene[2],  gene[3]とgene[4] ... のように親にする */

    /*** この部分を自分で書く  ***/

    c_pos = Random() * (G_LENGTH - 1) + 1;// Random() * (G_LENGTH - 1) + 1の計算結果をc_posにいれる。
    for (i = 1; i < POP_SIZE; i += 2) {//ループ処理前に1をいれ、POP_SIZEの処理を2ずつ繰り返す。
        for (j = c_pos; j < G_LENGTH; j++) {//ループ処理前に1をいれ、POP_SIZEの処理を2ずつ繰り返す。
            gene1[j] = gene[i][j]; //gene[i][j]をgene1[j]にコピーする。
            gene2[j] = gene[i + 1][j]; //gene[i+1][j]をgene2[j]にコピーする。
        }
        for (j = c_pos; j < G_LENGTH; j++) {//ループ処理前にc_posをいれ、G_LEMGTHの処理を繰り返す。
            gene[i][j] = gene2[j]; //gene2[j]をgene[i][j]にコピーする。
            gene[i + 1][j] = gene1[j]; //gene1[j]をgene[i+1][j]にコピーする。
        }
    }
}

/**********************************************************************
二点交叉 (余裕があれば)
引数 gene[p][i] : 遺伝子pのi番目の成分
**********************************************************************/
void two_crossover(int gene[POP_SIZE][G_LENGTH])
{
    int gene1[G_LENGTH]; /* 親1の遺伝子型 */
    int gene2[G_LENGTH]; /* 親2の遺伝子型 */
    int p, i;
    int c_pos1, c_pos2;  /* 交叉位置 (1 <= c_pos1,2 <= G_LENGTH-1) */
    int tmp;

    /* 交叉位置を1～G_LENGTH-1の範囲でランダムに2つ決め、その間を入れ替える。
    gene[1]とgene[2],  gene[3]とgene[4] ... のように親にする */

    /*** この部分を自分で書く  ***/
    c_pos1 = Random() * (G_LENGTH - 1) + 1; // Random() * (G_LENGTH - 1) + 1の計算結果をc_pos1にいれる。
    c_pos2 = Random() * (G_LENGTH - 1) + 1; //Random() * (G_LENGTH - 1) + 1の計算結果をc_pos2にいれる。
    if (c_pos1 > c_pos2) { tmp = c_pos1; c_pos2 = c_pos1; c_pos1 = tmp; } //もしc_pos1がc_pos2より大きければ、c_pos1をtmpにコピーして、c_pos1をc_pos2にコピーして、tmpをc_pos1にコピーする。

    for (p = 1; p < POP_SIZE; p += 2) {//ループ処理前に1をいれ、POP_SIZEの処理を2ずつ繰り返す。
        for (i = c_pos1; i < c_pos2; i++) {//ループ処理前にc_pos1をいれ、c_pos2の処理を繰り返す。
            gene1[i] = gene[p][i]; //gene[p][i]をgene1[i]にコピーする。
            gene2[i] = gene[p + 1][i];//gene[p+1][i]をgene2[i]にコピーする。
        }
        for (i = c_pos1; i < c_pos2; i++) {//ループ処理前にc_pos1をいれ、c_pos2の処理を繰り返す。
            gene[p][i] = gene2[i]; //gene[i]をgene[p][i]にコピーする。
            gene[p + 1][i] = gene1[i];//gene1[i]をgene[p+1][i]にコピーする。
        }
    }
}


/**********************************************************************
突然変異
引数 gene[p][i] : 遺伝子pのi番目の成分
**********************************************************************/
void mutation(int gene[POP_SIZE][G_LENGTH])
{
    int p, i;

    /* 0～1の乱数を発生させ、その値が M_RATE 以下ならば
    遺伝子の値をランダムに変える (0ならば1、1ならば0) */

    /*** この部分を自分で書く  ***/
    for (p = 1; p < POP_SIZE; p++)//ループ処理前に1をいれ、POP_SIZEの処理を繰り返す。
        for (i = 0; i < G_LENGTH; i++) //ループ処理前に0をいれ、G_LEMGTHの処理を繰り返す。
            if (Random() <= M_RATE) //もしRandom()はM_RATE以下なら、
                gene[p][i] = (gene[p][i] == 1) ? 0 : 1;//gene[p][i]が1なら0,1でなければ1。
}

/**********************************************************************
メインプログラム
**********************************************************************/
int main(int argc, char *argv[])
{
    int gene[POP_SIZE][G_LENGTH];
    double fitness[POP_SIZE];
    //int t;
    //FILE *fp;

    /* 適応度の変化を記録するファイルのオープン */
    if ((fp = fopen("result.dat", "w")) == NULL) {
        printf("Cannot open \"result.dat\"\n");
        exit(1);
    }

    /* シミュレーション条件の表示 */
    printf("個体数     : %d\n", POP_SIZE);
    printf("遺伝子長   : %d bit\n", G_LENGTH);
    printf("突然変異率 : %lf\n", M_RATE);


    init_gene(gene);              /* 遺伝子の初期化 */
    calc_fitness(gene, fitness);   /* 適応度の計算 */
    show_gene(t, gene, fitness, fp); /* 表示 */

    for (t = 1; t <= MAX_GEN; t++) {
        printf("<< 世代数 : %d >>\n", t);
        elite(gene, fitness);           /* エリート保存 */
        reproduction(gene, fitness);    /* ルーレット選択 */
        crossover(gene);               /* 単純交叉 */
        //two_crossover(gene);         /* 二点交叉 */
        mutation(gene);                /* 突然変異 */
        calc_fitness(gene, fitness);    /* 適応度の計算 */
        show_gene(t, gene, fitness, fp);  /* 表示 */
    }

    /* ファイルのクローズ */
    fclose(fp);

    return 0;
}
