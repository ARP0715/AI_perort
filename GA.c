//
// Created by Admin on 2018/07/25.
//

/********************************************************************
��`�I�A���S���Y��
********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <math.h>

#define POP_SIZE 5   /* �̐� (�K����ɐݒ�) */
#define G_LENGTH 10  /* �̂̈�`�q�^�̃r�b�g�� */
#define MAX_GEN  20  /* ���㐔 */
#define M_RATE   0.1 /* �ˑR�ψٗ� (0�`1) */

double avg_fit;  /* ���ϓK���x  */
double max_fit;
int t;
FILE *fp;


/********************************************************************
�����̔��� (Seed�̌���)
********************************************************************/
void init_rnd()
{
    srand((unsigned int)time(NULL));
}

/********************************************************************
�����̔��� (0�`1�̗���)
********************************************************************/
double Random()
{
    return((double)rand() / RAND_MAX);
}

/********************************************************************
��`�q�̏�����
���� gene[p][i] : ��`�qp��i�Ԗڂ̐���
********************************************************************/
void init_gene(int gene[POP_SIZE][G_LENGTH])
{
    int p, i;

    /* �����̏����� */
    init_rnd();

    /* ��`�q��������  0�`1�̗����𔭐����A0.5�ȏ�Ȃ�1
    0.5�����Ȃ�0 */
    printf("<< �����̌Q >>\n");

    /*** ���̕����������ŏ���  ***/
    for (p = 0; p < POP_SIZE; p++) {//���[�v�����O��0������APOP_SIZE�̏������J��Ԃ��B
        for (i = 0; i < G_LENGTH; i++) {//���[�v�����O��0������AG_LENGTH�̏������J��Ԃ��B

            if (Random() < 0.5) gene[p][i] = 0; //����Random()��0.5��菬�������0�ƂȂ�B
            else gene[p][i] = 1; //����ȊO��1�ƂȂ�B
        }
    }
}

/********************************************************************
�K���x�̌v�Z
���� gene[p][i] : ��`�qp��i�Ԗڂ̐���
fitness[p] : ��`�qp�̓K���x
********************************************************************/
void calc_fitness(int gene[POP_SIZE][G_LENGTH], double fitness[POP_SIZE])
{
    int p, i;

    /* �K���x�̌v�Z �O����5bit��0�̐� �㔼��5bit��1�̐� */

    /*** ���̕����������ŏ���  ***/
    for (p = 0; p < POP_SIZE; p++) {//���[�v�����O��0������APOP_SIZE�̏������J��Ԃ��B

        fitness[p] = 0;//fitness��0�������B
        for (i = 0; i < 5; i++) {//i��0����͂��܂�T�ɂȂ�܂ŌJ��Ԃ��B
            if (gene[p][i] == 0) fitness[p] += 1;// gene[p][i] == 0�Ȃ�fitness[p]��1�𑫂��B
        }
        for (; i < G_LENGTH; i++) {//i��G_LENGTH�̏������J��Ԃ��B
            if (gene[p][i] == 1) fitness[p] += 1;//����gene[p][i] == 1�Ȃ�fitness[p]��1�𑫂��B
        }
    }
}

/**********************************************************************
��`�q�̕\�� & �ő�K���x�E���ϓK���x�̌v�Z & �t�@�C���ւ̏����o��
���� t          : ���㐔
gene[p][i] : ��`�qp��i�Ԗڂ̐���
fitness[p] : ��`�qp�̓K���x
*fp        : �t�@�C���|�C���^
**********************************************************************/
void show_gene(int t, int gene[POP_SIZE][G_LENGTH], double fitness[POP_SIZE], FILE *fp) {
    //int t, int gene[POP_SIZE][G_LENGTH], double fitness[POP_SIZE], FILE *fp
    int p, i;
    double avg_fit;  /* ���ϓK���x  */
    double max_fit;  /* �ő�K���x  */

    /* �̂̒l�A�K���x�̕\�� */

    /*** ���̕����������ŏ���  ***/
    avg_fit = 0.0;//avg_fit��0.0�������B
    max_fit = 0.0;//max_fit��0.0�������B
    /* ���ρE�ő�K���x�̌v�Z */
    for (p = 0; p < POP_SIZE; p++) {//���[�v�����O��0������APOP_SIZE�̏������J��Ԃ��B
        if (max_fit < fitness[p]) max_fit = fitness[p];//����max_fit��fitness[p]��菬������΁Amax_fit��fitness[p]���R�s�[����B
        avg_fit += fitness[p]; // fitness[p]����avg_fit�𑫂��B
    }
    avg_fit /= POP_SIZE;//POP_SIZE����avg_fit������B



    /* ���ρE�ő�K���x�̕\�� */
    printf("���ϓK���x : %lf\n", avg_fit);
    printf("�ő�K���x : %lf\n", max_fit);

    /* ���ρE�ő�K���x���t�@�C���ɏ������� */
    fprintf(fp, "%d %lf %lf\n", t, avg_fit, max_fit);

}

/**********************************************************************
�̔ԍ� p1 �� p2 �̓K���x�ƈ�`�q������
���� p1, p2     : ��`�q�̔ԍ�
gene[p][i] : ��`�qp��i�Ԗڂ̐���
fitness[p] : ��`�qp�̓K���x
**********************************************************************/
void swap_gene(int p1, int p2, int gene[POP_SIZE][G_LENGTH],
               double fitness[POP_SIZE])
{
    int tmp;
    double f;
    int i;

    /* ��`�q�^�̌��� (��`�qp1�ƈ�`�qp2�̒l�����ւ���) */

    /*** ���̕����������ŏ���  ***/
    for (i = 0; i < G_LENGTH; i++) { //���[�v�����O��0������AG_LENGTH�̏������J��Ԃ��B
        tmp = gene[p1][i];//gene[p1][i]��tmp�ɃR�s�[����B
        gene[p1][i] = gene[p2][i];//gene[p2][i]��gene[p1][i]�ɃR�s�[����B
        gene[p2][i] = tmp; //tmp��gene[p2][i]�ɃR�s�[����B
    }

    /* �K���x�̌��� (��`�qp1�ƈ�`�qp2�̓K���x�̒l�����ւ���) */

    /*** ���̕����������ŏ���  ***/
    f = fitness[p1];//f��fitness[p1]�������B
    fitness[p1] = fitness[p2];//fitness[p2]��fitness[p1]�ɃR�s�[����B
    fitness[p2] = f;//f��fitness[p2]�ɃR�s�[����B


}

/**********************************************************************
�̔ԍ� p1 �̓K���x�ƈ�`�q�^�� p2 �ɃR�s�[
���� p1, p2     : ��`�q�̔ԍ�
gene[p][i] : ��`�qp��i�Ԗڂ̐���
fitness[p] : ��`�qp�̓K���x
**********************************************************************/
void copy_gene(int p1, int p2, int gene[POP_SIZE][G_LENGTH],
               double fitness[POP_SIZE])
{
    int i;

    /* ��`�q�̃R�s�[ (��`�qp1����`�qp2�ɃR�s�[����) */

    /*** ���̕����������ŏ���  ***/
    for (i = 0; i < G_LENGTH; i++) {//���[�v�����O��0������AG_LENGTH�̏������J��Ԃ��B
        gene[p2][i] = gene[p1][i]; //gene[p1][i]��gene[p2][i]�ɃR�s�[����B
    }

    /* �K���x�̃R�s�[ (��`�qp1�̓K���x����`�qp2�̓K���x�ɃR�s�[����)*/

    /*** ���̕����������ŏ���  ***/
    fitness[p2] = fitness[p1];//fitness[p1]��fitness[p2]�ɃR�s�[����B

}

/**********************************************************************
�G���[�g�ۑ�
(�ŏ��K���x�̌̂ɍő�K���x�̌̂̃f�[�^���R�s�[)
���� gene[p][i] : ��`�qp��i�Ԗڂ̐���
fitness[p] : ��`�qp�̓K���x
**********************************************************************/
void elite(int gene[POP_SIZE][G_LENGTH], double fitness[POP_SIZE])
{
    int p, i;
    double max_fitness = fitness[0];
    double min_fitness = fitness[0];
    int max_p = 0;
    int min_p = 0;

    /* �ő�K���x�̌�(max_p)�ƍŏ��K���x�̌�(min_p)�������� */

    /*** ���̕����������ŏ���  ***/
    for (p = 1; p < POP_SIZE; p++) {//���[�v�����O��1������APOP_SIZE�̏������J��Ԃ��B
        if (max_fitness < fitness[p]) {  //����max_fitness��fitness[p]��肿��������΁A
            max_fitness = fitness[p]; //fitness[p]��max_fitness�ɃR�s�[����B
            max_p = p;//p��max_p�ɃR�s�[����B
        }
        if (min_fitness > fitness[p]) {//����min_fitness��fitness[p]���傫����΁A
            min_fitness = fitness[p];//fitness[p]��min_fitness�ɃR�s�[����B
            min_p = p;//p��min_p�ɃR�s�[����B
        }
    }



    /* �ŏ��K���x�̌̂ɍő�K���x�̌̂��R�s�[ */
    copy_gene(max_p, min_p, gene, fitness);
    /* �ő�K���x�̌̂�0�ԖڂɈړ� */
    swap_gene(0, max_p, gene, fitness);
}

/**********************************************************************
���[���b�g�I��
���� gene[p][i] : ��`�qp��i�Ԗڂ̐���
fitness[p] : ��`�qp�̓K���x
**********************************************************************/
void reproduction(int gene[POP_SIZE][G_LENGTH], double fitness[POP_SIZE])
{
    double sum_of_fitness; /* �̂̓K���x�̑��a */
    double border;         /* ���[���b�g��̌̊Ԃ̋��E */
    double r;              /* ���[���b�g��̑I���ʒu */
    int p, i;               /* �I�΂ꂽ�̂̔ԍ� */
    int num;               /* 0 <= num <= POP_SIZE-1 */
    int new_gene[POP_SIZE][G_LENGTH];

    /* ���[���b�g��1���� sum_of_fitness �����߂� */
    sum_of_fitness = 0.0;
    for (p = 0; p<POP_SIZE; p++) {
        sum_of_fitness += fitness[p];
    }

    /* ���[���b�g�� POP_SIZE �񂾂��񂵂Ď�����̌̂�I�� */
    for (p = 1; p<POP_SIZE; p++) {
        /* ���[���b�g���񂵂ďꏊ��I��
        r : �I�΂ꂽ�ʒu (0 <= r <= sum_of_fitness) */
        r = sum_of_fitness*Random();
        /* �I�΂ꂽ�ꏊ�ɊY������̂����Ԃ����ׂ�
        num : �I�΂ꂽ�̂̔ԍ� (0 <= num <= POP_SIZE-1) */
        num = 0;
        border = fitness[0]; /* �̊Ԃ̋��E */
        while (border<r) {
            num++;
            border += fitness[num];
        }

        /* ��`�q�̑�� */
        for (i = 0; i<G_LENGTH; i++) {
            new_gene[p][i] = gene[num][i];
        }
    }

    /* ��`�q�̃R�s�[ */
    for (p = 1; p<POP_SIZE; p++) {
        for (i = 0; i<G_LENGTH; i++) {
            gene[p][i] = new_gene[p][i];
        }
    }
}

/**********************************************************************
��_����
���� gene[p][i] : ��`�qp��i�Ԗڂ̐���
**********************************************************************/
void crossover(int gene[POP_SIZE][G_LENGTH])
{
    int gene1[G_LENGTH]; /* �e1�̈�`�q�^ */
    int gene2[G_LENGTH]; /* �e2�̈�`�q�^ */
    int i, j;
    int c_pos;   /* �����ʒu (1 <= c_pos <= G_LENGTH-1) */

    /* �����ʒu��1�`G_LENGTH-1�͈̔͂Ń����_���Ɍ��߁A
    ������������ւ���B
    gene[1]��gene[2],  gene[3]��gene[4] ... �̂悤�ɐe�ɂ��� */

    /*** ���̕����������ŏ���  ***/

    c_pos = Random() * (G_LENGTH - 1) + 1;// Random() * (G_LENGTH - 1) + 1�̌v�Z���ʂ�c_pos�ɂ����B
    for (i = 1; i < POP_SIZE; i += 2) {//���[�v�����O��1������APOP_SIZE�̏�����2���J��Ԃ��B
        for (j = c_pos; j < G_LENGTH; j++) {//���[�v�����O��1������APOP_SIZE�̏�����2���J��Ԃ��B
            gene1[j] = gene[i][j]; //gene[i][j]��gene1[j]�ɃR�s�[����B
            gene2[j] = gene[i + 1][j]; //gene[i+1][j]��gene2[j]�ɃR�s�[����B
        }
        for (j = c_pos; j < G_LENGTH; j++) {//���[�v�����O��c_pos������AG_LEMGTH�̏������J��Ԃ��B
            gene[i][j] = gene2[j]; //gene2[j]��gene[i][j]�ɃR�s�[����B
            gene[i + 1][j] = gene1[j]; //gene1[j]��gene[i+1][j]�ɃR�s�[����B
        }
    }
}

/**********************************************************************
��_���� (�]�T�������)
���� gene[p][i] : ��`�qp��i�Ԗڂ̐���
**********************************************************************/
void two_crossover(int gene[POP_SIZE][G_LENGTH])
{
    int gene1[G_LENGTH]; /* �e1�̈�`�q�^ */
    int gene2[G_LENGTH]; /* �e2�̈�`�q�^ */
    int p, i;
    int c_pos1, c_pos2;  /* �����ʒu (1 <= c_pos1,2 <= G_LENGTH-1) */
    int tmp;

    /* �����ʒu��1�`G_LENGTH-1�͈̔͂Ń����_����2���߁A���̊Ԃ����ւ���B
    gene[1]��gene[2],  gene[3]��gene[4] ... �̂悤�ɐe�ɂ��� */

    /*** ���̕����������ŏ���  ***/
    c_pos1 = Random() * (G_LENGTH - 1) + 1; // Random() * (G_LENGTH - 1) + 1�̌v�Z���ʂ�c_pos1�ɂ����B
    c_pos2 = Random() * (G_LENGTH - 1) + 1; //Random() * (G_LENGTH - 1) + 1�̌v�Z���ʂ�c_pos2�ɂ����B
    if (c_pos1 > c_pos2) { tmp = c_pos1; c_pos2 = c_pos1; c_pos1 = tmp; } //����c_pos1��c_pos2���傫����΁Ac_pos1��tmp�ɃR�s�[���āAc_pos1��c_pos2�ɃR�s�[���āAtmp��c_pos1�ɃR�s�[����B

    for (p = 1; p < POP_SIZE; p += 2) {//���[�v�����O��1������APOP_SIZE�̏�����2���J��Ԃ��B
        for (i = c_pos1; i < c_pos2; i++) {//���[�v�����O��c_pos1������Ac_pos2�̏������J��Ԃ��B
            gene1[i] = gene[p][i]; //gene[p][i]��gene1[i]�ɃR�s�[����B
            gene2[i] = gene[p + 1][i];//gene[p+1][i]��gene2[i]�ɃR�s�[����B
        }
        for (i = c_pos1; i < c_pos2; i++) {//���[�v�����O��c_pos1������Ac_pos2�̏������J��Ԃ��B
            gene[p][i] = gene2[i]; //gene[i]��gene[p][i]�ɃR�s�[����B
            gene[p + 1][i] = gene1[i];//gene1[i]��gene[p+1][i]�ɃR�s�[����B
        }
    }
}


/**********************************************************************
�ˑR�ψ�
���� gene[p][i] : ��`�qp��i�Ԗڂ̐���
**********************************************************************/
void mutation(int gene[POP_SIZE][G_LENGTH])
{
    int p, i;

    /* 0�`1�̗����𔭐������A���̒l�� M_RATE �ȉ��Ȃ��
    ��`�q�̒l�������_���ɕς��� (0�Ȃ��1�A1�Ȃ��0) */

    /*** ���̕����������ŏ���  ***/
    for (p = 1; p < POP_SIZE; p++)//���[�v�����O��1������APOP_SIZE�̏������J��Ԃ��B
        for (i = 0; i < G_LENGTH; i++) //���[�v�����O��0������AG_LEMGTH�̏������J��Ԃ��B
            if (Random() <= M_RATE) //����Random()��M_RATE�ȉ��Ȃ�A
                gene[p][i] = (gene[p][i] == 1) ? 0 : 1;//gene[p][i]��1�Ȃ�0,1�łȂ����1�B
}

/**********************************************************************
���C���v���O����
**********************************************************************/
int main(int argc, char *argv[])
{
    int gene[POP_SIZE][G_LENGTH];
    double fitness[POP_SIZE];
    //int t;
    //FILE *fp;

    /* �K���x�̕ω����L�^����t�@�C���̃I�[�v�� */
    if ((fp = fopen("result.dat", "w")) == NULL) {
        printf("Cannot open \"result.dat\"\n");
        exit(1);
    }

    /* �V�~�����[�V���������̕\�� */
    printf("�̐�     : %d\n", POP_SIZE);
    printf("��`�q��   : %d bit\n", G_LENGTH);
    printf("�ˑR�ψٗ� : %lf\n", M_RATE);


    init_gene(gene);              /* ��`�q�̏����� */
    calc_fitness(gene, fitness);   /* �K���x�̌v�Z */
    show_gene(t, gene, fitness, fp); /* �\�� */

    for (t = 1; t <= MAX_GEN; t++) {
        printf("<< ���㐔 : %d >>\n", t);
        elite(gene, fitness);           /* �G���[�g�ۑ� */
        reproduction(gene, fitness);    /* ���[���b�g�I�� */
        crossover(gene);               /* �P������ */
        //two_crossover(gene);         /* ��_���� */
        mutation(gene);                /* �ˑR�ψ� */
        calc_fitness(gene, fitness);    /* �K���x�̌v�Z */
        show_gene(t, gene, fitness, fp);  /* �\�� */
    }

    /* �t�@�C���̃N���[�Y */
    fclose(fp);

    return 0;
}
