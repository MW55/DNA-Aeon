pipeline {
    //agent any
    agent {
        docker {
                image 'ubuntu:20.04' // 'qnib/pytest'
        }
    }

    environment {
        CODE_COV_TOKEN = credentials('CODE_COV_TOKEN')
    }

    stages {
        stage('Build') {
            steps {
                sh script: '''DEBIAN_FRONTEND=noninteractive apt update &&DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata && DEBIAN_FRONTEND=noninteractive apt install -y python3.8 apt-utils python3-pip lsb-release wget curl software-properties-common'''
                sh script: '''bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)" && python3.8 -m pip install --upgrade pip && python3.8 -m pip install pytest && python3.8 -m pip install pytest-cov && \
                 python3.8 -m pip install codecov && python3.8 -m pip install numpy  && python3.8 -m pip install pytest-cov && python3.8 -m pip install -r requirements.txt && python3.8 setup.py install'''
            }
        }
        stage('Test') {
            steps {
                //withCredentials([file(credentialsId: 'secret', variable: 'CODE_COV_TOKEN')]) {
                withCredentials([string(credentialsId: 'CODE_COV_TOKEN', variable: 'CODE_COV_TOKEN')]) {
                    sh script: '''
                    python3.8 -m pytest tests/ --cov=./norec4dna && \
                      curl -s https://codecov.io/bash | bash -s - -t $CODE_COV_TOKEN'''
                }
            }
        }
        stage('Deploy') {
            steps {
                echo 'Deploying....'
            }
        }
    }
}